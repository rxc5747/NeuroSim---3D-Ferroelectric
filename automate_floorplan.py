import os
import re
import csv
import subprocess
import sys
import signal
import time
 
# --- CONFIGURATION ---
model_name = "resnet18"
memcelltype = 5                                   # 5 = 2TnC, 6 = 1TnC
n_values = [4, 8, 16, 32, 64]
subarray_shapes = [(128, 128), (256, 256), (512, 512)]   # (physical rows, cols)
integration_modes = ["CNA", "CUA", "CBA"]          # typedef.h enum names
# ---------------------
 
def modify_file(filepath, pattern, replacement):
    """Utility to find and replace text in a file."""
    with open(filepath, "r") as f:
        content = f.read()
 
    # Replace the target string (warn if the pattern is missing, so a
    # silent no-op regex can't corrupt the whole sweep)
    new_content, count = re.subn(pattern, replacement, content)
    if count == 0:
        print(f"   [WARN] pattern not found in {filepath}: {pattern}")
 
    with open(filepath, "w") as f:
        f.write(new_content)
 
 
 
def parse_area_lines(fname):
    import re as _re
    F = r"([-+]?[\d.]+(?:[eE][-+]?\d+)?)"
    keys = [("chip_area_um2",      r"ChipArea\s*:?\s*" + F),
            ("staircase_um2",      r"\[3D\] staircase / via routing:\s*" + F),
            ("memory_die_um2",     r"\[3D\] memory / logic die:\s*" + F),
            ("logic_die_um2",      r"\[3D\] memory / logic die:\s*[-+\d.eE]+\s*/\s*" + F),
            ("bond_pad_um2",       r"\[3D\] bond pad / chip:\s*" + F),
            ("chip_footprint_um2", r"\[3D\] bond pad / chip:\s*[-+\d.eE]+\s*/\s*" + F),
            ("array_eff_pct",      r"\[3D\] array efficiency:\s*" + F)]
    vals = {k: "" for k, _ in keys}
    try:
        text = open(fname, errors="ignore").read()
        for k, pat in keys:
            m = _re.search(pat, text)
            if m: vals[k] = m.group(1)
    except FileNotFoundError:
        pass
    return vals
 
def run_is_good(fname):
    """A result file counts as successful if it exists and reached the
    chip-level summary (failed runs only contain the ERROR + segfault)."""
    try:
        return "ChipArea" in open(fname, errors="ignore").read()
    except FileNotFoundError:
        return False
 
summary_rows = []
 
for n in n_values:
    for (rows, cols) in subarray_shapes:
        expected = [f"results_{model_name}_2T{n}C_{rows}x{cols}_{m}_no_pipeline.txt" for m in integration_modes]
        if all(run_is_good(f) for f in expected):
            print(f"        SKIPPING 2T{n}C {rows}x{cols} -- all 3 modes already complete")
            for m, f in zip(integration_modes, expected):
                row = {"n": n, "rows": rows, "cols": cols, "mode": m, "file": f}
                row.update(parse_area_lines(f))
                summary_rows.append(row)
            continue
 
        print(f"        STARTING PIPELINE FOR 2T{n}C  {rows}x{cols}              ")
 
        # 1. Calculate the physics bounds
        logical_rows = n * rows
        pe_size = logical_rows * 2
        tile_size = logical_rows * 4
 
        # 2. Modify Param.cpp
        print("-> Updating Param.cpp...")
        modify_file("NeuroSIM/Param.cpp", r"memcelltype\s*=\s*\d+;", f"memcelltype = {memcelltype};")
        modify_file("NeuroSIM/Param.cpp", r"bitsPerCell\s*=\s*\d+;", f"bitsPerCell = {n};")
        modify_file("NeuroSIM/Param.cpp", r"numRowSubArrayPhysical\s*=\s*\d+;", f"numRowSubArrayPhysical = {rows};")
        modify_file("NeuroSIM/Param.cpp", r"numColSubArray\s*=\s*\d+;", f"numColSubArray = {cols};")
 
        # 3. Modify main.cpp -- REQUIRED above ~2048 logical rows, otherwise
        #    ChipDesignInitialize fails with "SubArray Size is too large"
        print("-> Updating main.cpp...")
        modify_file("NeuroSIM/main.cpp", r"maxPESizeNM\s*=\s*\d+;", f"maxPESizeNM = {pe_size};")
        modify_file("NeuroSIM/main.cpp", r"maxTileSizeCM\s*=\s*\d+;", f"maxTileSizeCM = {tile_size};")
 
        # 4. Modify inference.py  (ACTIVE now: the shape sweep must reach the
        #    Python-side mapping; generalized from your [\d+, 512] pattern)
        print("-> Updating inference.py...")
        modify_file("inference.py", r"args\.sub_array\s*=\s*\[\d+,\s*\d+\]", f"args.sub_array = [{logical_rows}, {cols}]")
 
        # 5. Compile NeuroSim (C++)
        print("-> Compiling NeuroSim (make clean && make)...")
        subprocess.run("make clean", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
        subprocess.run("make", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
 
        # 6. Run Python Inference with Auto-Interrupt
        #    (one inference per (n, shape); integrationMode is C++-only, so the
        #    traces are reused across CNA/CUA/CBA below)
        print("-> Running inference.py (Waiting for 4% batch completion)...")
        inf_cmd = f"python inference.py --dataset cifar10 --model {model_name} --data_path ./datasets --mem_type \"capacitive\""
        # NOTE: kept "capacitive" to match your validated pipeline; for the 2TnC
        # current-domain read you may want "resistive" -- change in one place here.
 
        # Start the process in a new process group so we can cleanly kill it later
        process = subprocess.Popen(inf_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, preexec_fn=os.setsid)
 
        hardware_sim_started = False
 
        # Listen to the terminal output line-by-line
        for line in process.stdout:
            print(line, end="")  # Print to terminal so you can monitor it
 
            # Step 1: Detect that the final simulation phase has begun
            if "Running inference with detailed hardware simulation" in line:
                hardware_sim_started = True
 
            # Step 2: Once started, wait for the first batch to finish (4% or 1/25)
            if hardware_sim_started and ("4%|" in line or "1/25" in line or "800 test images" in line):
                print("\n-> [Auto-Interrupt] 4% batch mark reached! CSVs should be fully saved.")
                time.sleep(2)  # Brief buffer to ensure file I/O is fully written to disk
 
                print("-> [Auto-Interrupt] Sending Ctrl+C (SIGTERM)...")
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)  # Safely kill the process tree
                break
 
        process.wait()  # Wait for Python to close
 
        # 7. Sweep the three area configurations on the SAME traces:
        #    only Param.cpp changes -> recompile -> re-run the trace command
        for mode in integration_modes:
            if run_is_good(f"results_{model_name}_2T{n}C_{rows}x{cols}_{mode}_no_pipeline.txt"):
                print(f"-> integrationMode = {mode}: already complete, skipping")
                continue
            print(f"-> integrationMode = {mode}: recompiling...")
            modify_file("NeuroSIM/Param.cpp", r"integrationMode\s*=\s*\w+;", f"integrationMode = {mode};")
            subprocess.run("make clean", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
            subprocess.run("make", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
 
            print("-> Running NeuroSim Hardware Trace...")
            log_filename = f"results_{model_name}_2T{n}C_{rows}x{cols}_{mode}_no_pipeline.txt"
            trace_cmd = f"bash ./layer_record_{model_name}/trace_command.sh"
 
            with open(log_filename, "w") as log_file:
                subprocess.run(trace_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
 
            print(f"-> DONE! Saved hardware performance to {log_filename}\n")
 
            # 8. parse the area lines into the running CSV summary
            row = {"n": n, "rows": rows, "cols": cols, "mode": mode, "file": log_filename}
            row.update(parse_area_lines(log_filename))
            summary_rows.append(row)
            with open("sweep_summary.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
                w.writeheader()
                w.writerows(summary_rows)
 
print("ALL RUNS COMPLETED SUCCESSFULLY.")
