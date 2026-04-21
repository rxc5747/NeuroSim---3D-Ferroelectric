import os
import re
import subprocess
import sys
import signal
import time

# --- CONFIGURATION ---
model_name = "resnet18"
n_values = [4, 8, 16, 32, 64]
# ---------------------

def modify_file(filepath, pattern, replacement):
    """Utility to find and replace text in a file."""
    with open(filepath, "r") as f:
        content = f.read()
    
    # Replace the target string
    new_content = re.sub(pattern, replacement, content)
    
    with open(filepath, "w") as f:
        f.write(new_content)

for n in n_values:
    print(f"==================================================")
    print(f"        STARTING PIPELINE FOR 1T{n}C              ")
    print(f"==================================================")
    
    # 1. Calculate the physics bounds
    logical_rows = n * 512
    pe_size = logical_rows * 2
    tile_size = logical_rows * 4
    
    # 2. Modify Param.cpp
    print("-> Updating Param.cpp...")
    modify_file("NeuroSIM/Param.cpp", r"bitsPerCell\s*=\s*\d+;", f"bitsPerCell = {n};")
    
    # # 3. Modify main.cpp
    # print("-> Updating main.cpp...")
    # modify_file("NeuroSIM/main.cpp", r"maxPESizeNM\s*=\s*\d+;", f"maxPESizeNM = {pe_size};")
    # modify_file("NeuroSIM/main.cpp", r"maxTileSizeCM\s*=\s*\d+;", f"maxTileSizeCM = {tile_size};")
    # 
    # # 4. Modify inference.py
    # print("-> Updating inference.py...")
    # modify_file("inference.py", r"args\.sub_array\s*=\s*\[\d+,\s*512\]", f"args.sub_array = [{logical_rows}, 512]")
    
    # 5. Compile NeuroSim (C++)
    print("-> Compiling NeuroSim (make clean && make)...")
    subprocess.run("make clean", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
    subprocess.run("make", shell=True, cwd="NeuroSIM", stdout=subprocess.DEVNULL)
    
    # 6. Run Python Inference with Auto-Interrupt
    print("-> Running inference.py (Waiting for 4% batch completion)...")
    inf_cmd = f"python inference.py --dataset cifar10 --model {model_name} --data_path ./datasets --mem_type \"capacitive\""
    
    # Start the process in a new process group so we can cleanly kill it later
    process = subprocess.Popen(inf_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, preexec_fn=os.setsid)
    
    hardware_sim_started = False

    # Listen to the terminal output line-by-line
    for line in process.stdout:
        print(line, end="") # Print to terminal so you can monitor it
        
        # Step 1: Detect that the final simulation phase has begun
        if "Running inference with detailed hardware simulation" in line:
            hardware_sim_started = True
            
        # Step 2: Once started, wait for the first batch to finish (4% or 1/25)
        if hardware_sim_started and ("4%|" in line or "1/25" in line or "800 test images" in line):
            print("\n-> [Auto-Interrupt] 4% batch mark reached! CSVs should be fully saved.")
            time.sleep(2) # Brief buffer to ensure file I/O is fully written to disk
            
            print("-> [Auto-Interrupt] Sending Ctrl+C (SIGTERM)...")
            os.killpg(os.getpgid(process.pid), signal.SIGTERM) # Safely kill the process tree
            break
            
    process.wait() # Wait for Python to close
    
    # 7. Execute the bash trace command and save output to a log file
    print("-> Running NeuroSim Hardware Trace...")
    log_filename = f"results_{model_name}_1T{n}C_no_pipeline.txt" 
    trace_cmd = f"bash ./layer_record_{model_name}/trace_command.sh"
    
    with open(log_filename, "w") as log_file:
        subprocess.run(trace_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        
    print(f"-> DONE! Saved hardware performance to {log_filename}\n")

print("ALL RUNS COMPLETED SUCCESSFULLY.")
