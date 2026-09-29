#!/usr/bin/env python3
"""
automate.py -- full CIM sweep over the config grid

    n         in {64, 256}                      (stacked capacitor planes = bitsPerCell)
    rows x cols in {64x64, 64x512, 512x64, 512x512}

    => 8 runs.

WHAT CHANGED vs THE ORIGINAL SCRIPT
  1. Every source edit is VERIFIED.  re.sub() silently does nothing when the
     pattern misses, so the original script would happily run all 8 configs
     against an unmodified Param.cpp and hand you 8 identical result files with
     no error.  This version aborts the moment a pattern fails to match, and
     prints the old -> new value for every edit it makes.
  2. The `make` return code is CHECKED.  The original ignored it, so a
     compile error left the previous binary in place and the run silently
     reported the PREVIOUS config's numbers.  That is the worst failure mode
     here because the output looks completely plausible.
  3. Sources are backed up before the sweep and restored afterwards (and on
     Ctrl-C / crash), so the repo never ends up stuck on the last config.
  4. The auto-interrupt has a watchdog timeout.  The original loop blocks
     forever if the "4%" marker never appears.
  5. Logs are per-config and the script skips a config whose log already
     exists, so an interrupted sweep resumes instead of restarting.
  6. A summary table at the end tells you which runs produced output and which
     did not.

ROW SEMANTICS -- READ THIS ONCE
  For _2TnC / _1TnC, NeuroSim carries two row counts:
     numRowSubArrayPhysical = physical pillar rows          (subArray->numRow)
     numRowSubArray         = physical * bitsPerCell        (logical / mapping)
  The grid in your spreadsheet is read as PHYSICAL rows x columns, matching
  your original script (which used 512 physical rows and 512 columns and set
  logical_rows = n * 512).  If your table actually means logical rows, flip
  ROWS_ARE_PHYSICAL to False and the derivation inverts.
  The script prints the derived values for every config, so check the first one
  before letting it run for hours.

USAGE
    python3 automate.py                 # run the sweep
    python3 automate.py --dry-run       # apply+verify the edits only, no build/run
    python3 automate.py --only 64x512   # run a subset (repeatable)
"""

import argparse
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
MODEL_NAME = "resnet18"
CELL_LABEL = "2T{n}C"                 # used in the log file name

N_VALUES   = [64, 256]
GEOMETRIES = [(64, 64), (64, 512), (512, 64), (512, 512)]   # (rows, cols)

ROWS_ARE_PHYSICAL = True              # the grid's "rows" are physical pillar rows

# What goes into Param.cpp numRowSubArray -- the number the MAPPER sees.
#   "physical" : numRowSubArray = physical rows.  n lives in bitsPerCell and
#                therefore in the array/staircase/plane-decoder models only.
#                This is the only setting that RUNS for ResNet-18 (see the
#                preflight below) and it matches the memory-mode access
#                definition (one access = one capacitor plane).
#   "logical"  : numRowSubArray = physical * n.  The mapper then believes the
#                subarray holds n x more weight rows.  NeuroSim's chip
#                hierarchy cannot express this for ResNet-18 -- every config
#                trips the ChipFloorPlan guard and segfaults.
MAPPER_ROWS = "physical"

# parallelRead (argv[5] of ./NeuroSIM/main, and --parallel_read in the python).
# None -> leave whatever inference.py already has.  A number -> force it.
# Keep it FIXED across the sweep, or the ADC levelOutput (and therefore the ADC
# peripheral cost) changes between configs and the comparison is not controlled.
PARALLEL_READ = 64

PARAM_CPP     = "NeuroSIM/Param.cpp"
MAIN_CPP      = "NeuroSIM/main.cpp"
INFERENCE_PY  = "inference.py"

INFER_CMD = ("python inference.py --dataset cifar10 --model {model} "
             "--data_path ./datasets --mem_type \"capacitive\"")

# how long to wait for the hardware-sim marker before giving up on a config
INFER_TIMEOUT_S = 60 * 60 * 6          # 6 h
# extra seconds to let CSV writes flush after the 4% mark
FLUSH_WAIT_S = 3

LOG_DIR = "sweep_logs"

# Files the script rewrites.  Backed up before the sweep, restored after.
TOUCHED_FILES = [PARAM_CPP, MAIN_CPP, INFERENCE_PY]


# --------------------------------------------------------------------------
# edit helpers -- these are the whole point of this rewrite
# --------------------------------------------------------------------------
class EditError(RuntimeError):
    pass


def edit(filepath, pattern, replacement, label, required=True):
    """Regex-replace in a file and PROVE it happened.

    Returns True if a substitution was made.  Raises EditError if `required`
    and the pattern matched nothing -- which is the failure the original
    script swallowed.
    """
    if not os.path.exists(filepath):
        if required:
            raise EditError(f"{label}: {filepath} does not exist")
        print(f"     skip  {label}: {filepath} not found")
        return False

    with open(filepath, "r") as f:
        content = f.read()

    matches = re.findall(pattern, content)
    if not matches:
        if required:
            raise EditError(
                f"{label}: pattern did not match anything in {filepath}\n"
                f"         pattern: {pattern}\n"
                f"         Fix the pattern (or the file) -- do NOT let the sweep\n"
                f"         continue, every run would use the same stale value."
            )
        print(f"     skip  {label}: no match in {filepath}")
        return False

    new_content, count = re.subn(pattern, replacement, content)
    if new_content == content:
        print(f"     ok    {label}: already set ({count} site(s))")
        return True

    with open(filepath, "w") as f:
        f.write(new_content)
    print(f"     set   {label}: {count} site(s) -> {replacement}")
    return True


def read_back(filepath, pattern, label):
    """Read a value straight out of the file after editing, to confirm."""
    with open(filepath, "r") as f:
        m = re.search(pattern, f.read())
    return m.group(1) if m else "<not found>"


# --------------------------------------------------------------------------
# backup / restore
# --------------------------------------------------------------------------
def backup_sources():
    saved = {}
    for p in TOUCHED_FILES:
        if os.path.exists(p):
            b = p + ".sweep_orig"
            shutil.copyfile(p, b)
            saved[p] = b
    print(f"-> backed up {len(saved)} source file(s)")
    return saved


def restore_sources(saved):
    for p, b in saved.items():
        if os.path.exists(b):
            shutil.copyfile(b, p)
            os.remove(b)
    if saved:
        print(f"-> restored {len(saved)} source file(s) to their original state")


# --------------------------------------------------------------------------
# chip-hierarchy preflight -- mirrors Chip.cpp exactly
# --------------------------------------------------------------------------
def _net_path():
    for p in (f"NeuroSIM/NetWork_{MODEL_NAME}.csv", f"NetWork_{MODEL_NAME}.csv"):
        if os.path.exists(p):
            return p
    return None


def hierarchy_limits(numRowSubArray, synapseBit=8, cellBit=1, numRowPerSynapse=1):
    """Recompute maxPESizeNM / maxTileSizeCM the way ChipDesignInitialize does.

    Chip.cpp:130/135   minCube = 2^ceil(log2(OFM_ch * numColPerSynapse))
                       maxPESizeNM / maxTileSizeCM = max over the marked layers
    Chip.cpp:214       novel mapping  : ERROR if maxPESizeNM  < 2 * numRowSubArray
    Chip.cpp:274       conventional   : ERROR if maxTileSizeCM < 4 * numRowSubArray
    After the ERROR the desired* outputs stay 0 and the caller dereferences
    them -> segmentation fault.  That is the core dump, not an OOM.
    """
    p = _net_path()
    if p is None:
        return None
    rows = []
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            v = [float(x) for x in line.split(",")]
            rows.append(v)
    numColPerSynapse = math.ceil(synapseBit / cellBit)
    cube = lambda ofm: 2 ** math.ceil(math.log2(ofm * numColPerSynapse))

    kk = [int(r[3] * r[4]) for r in rows if int(r[3] * r[4]) != 1]
    numPENM = max(set(kk), key=kk.count) if kk else 0

    peNM, tileCM = 0.0, 0.0
    for r in rows:
        ifm, k1, k2, ofm = r[2], r[3], r[4], r[5]
        isNM = (int(k1 * k2) == numPENM) and (ifm * k1 * k2 * numRowPerSynapse >= numRowSubArray)
        if isNM:
            peNM = max(peNM, cube(ofm))
        else:
            tileCM = max(tileCM, cube(ofm))
    return dict(maxPESizeNM=peNM, maxTileSizeCM=tileCM,
                nm_ok=not (peNM < 2 * numRowSubArray),
                cm_ok=not (tileCM < 4 * numRowSubArray),
                need_nm=2 * numRowSubArray, need_cm=4 * numRowSubArray)


def preflight(numRowSubArray):
    h = hierarchy_limits(numRowSubArray)
    if h is None:
        print("   preflight: NetWork csv not found, skipping hierarchy check")
        return True
    ok = h["nm_ok"] and h["cm_ok"]
    print(f"   preflight: maxPESizeNM={h['maxPESizeNM']:.0f} (needs >= {h['need_nm']:.0f})  "
          f"maxTileSizeCM={h['maxTileSizeCM']:.0f} (needs >= {h['need_cm']:.0f})  "
          f"-> {'OK' if ok else 'WILL SEGFAULT'}")
    if not ok:
        raise EditError(
            "chip hierarchy would reject numRowSubArray=%d.\n"
            "         Chip.cpp ChipFloorPlan requires\n"
            "             maxPESizeNM  >= 2 * numRowSubArray   (novel mapping)\n"
            "             maxTileSizeCM >= 4 * numRowSubArray  (conventional)\n"
            "         and both are derived from the WIDEST LAYER of the network,\n"
            "         not from anything you can set:\n"
            "             minCube = 2^ceil(log2(OFM_ch * numColPerSynapse))\n"
            "         For %s at 8-bit weights / 1 bit-per-cell that caps\n"
            "         numRowSubArray at 2048.  Above it the guard fires, the\n"
            "         desired* outputs stay 0, and main segfaults.\n"
            "         Fix: set MAPPER_ROWS = \"physical\" at the top of this file."
            % (numRowSubArray, MODEL_NAME))
    return ok


# --------------------------------------------------------------------------
# per-config pipeline
# --------------------------------------------------------------------------
def apply_config(n, rows, cols):
    """Write this config into the sources.  Raises EditError on any miss."""
    if ROWS_ARE_PHYSICAL:
        phys_rows    = rows
        logical_rows = rows * n
    else:
        logical_rows = rows
        phys_rows    = rows // n
        if phys_rows * n != logical_rows:
            raise EditError(f"logical rows {logical_rows} is not divisible by n={n}")

    pe_size   = logical_rows * 2
    tile_size = logical_rows * 4

    print(f"   physical rows {phys_rows} | logical rows {logical_rows} "
          f"| cols {cols} | bitsPerCell {n}")

    # ---- Param.cpp ----------------------------------------------------
    edit(PARAM_CPP, r"bitsPerCell\s*=\s*\d+;",
         f"bitsPerCell = {n};", "Param.bitsPerCell")

    mapper_rows = logical_rows if MAPPER_ROWS == "logical" else phys_rows
    print(f"   mapper sees numRowSubArray = {mapper_rows} "
          f"(MAPPER_ROWS = {MAPPER_ROWS!r})")
    preflight(mapper_rows)

    edit(PARAM_CPP, r"numRowSubArray\s*=\s*\d+;",
         f"numRowSubArray = {mapper_rows};", "Param.numRowSubArray")

    edit(PARAM_CPP, r"numColSubArray\s*=\s*\d+;",
         f"numColSubArray = {cols};", "Param.numColSubArray")

    # present only in the 3D FeRAM fork; not required so the script still
    # works on a stock tree
    edit(PARAM_CPP, r"numRowSubArrayPhysical\s*=\s*\d+;",
         f"numRowSubArrayPhysical = {phys_rows};",
         "Param.numRowSubArrayPhysical", required=False)

    # ---- main.cpp -----------------------------------------------------
    # DISABLED, and it matters why: maxPESizeNM and maxTileSizeCM are OUTPUT
    # parameters of ChipDesignInitialize() (Chip.cpp:98-99 zero them, :130/:135
    # fill them from the network's widest layer).  Patching a literal in
    # main.cpp is overwritten at runtime and does nothing.  Leaving the edits in
    # made the sweep look like it was controlling the hierarchy when it was not.
    # The hierarchy limit is enforced by the preflight check instead.

    # ---- inference.py -------------------------------------------------
    # The python mapper and the C++ model must agree on the subarray shape or
    # the whole run is meaningless.  Try the known spellings; require one.
    done = edit(INFERENCE_PY,
                r"args\.sub_array\s*=\s*\[\s*\d+\s*,\s*\d+\s*\]",
                f"args.sub_array = [{logical_rows}, {cols}]",
                "inference.args.sub_array", required=False)
    if not done:
        done = edit(INFERENCE_PY,
                    r"args\.subArray\s*=\s*\d+",
                    f"args.subArray = {logical_rows}",
                    "inference.args.subArray", required=False)
    if not done:
        raise EditError(
            "inference.py: could not find a sub-array assignment to patch.\n"
            "         Looked for 'args.sub_array = [R, C]' and 'args.subArray = N'.\n"
            "         Open inference.py, find the line that sets the subarray\n"
            "         shape, and add its pattern to apply_config().\n"
            "         Without this the python mapper and the C++ model disagree\n"
            "         and every number the sweep produces is meaningless."
        )

    if PARALLEL_READ is not None:
        edit(INFERENCE_PY, r"args\.parallel_read\s*=\s*\d+",
             f"args.parallel_read = {PARALLEL_READ}",
             "inference.args.parallel_read", required=False)
        edit(INFERENCE_PY, r"args\.parallelRead\s*=\s*\d+",
             f"args.parallelRead = {PARALLEL_READ}",
             "inference.args.parallelRead", required=False)

    # ---- read back, so the printed value is what is ON DISK -----------
    v_bits = read_back(PARAM_CPP, r"bitsPerCell\s*=\s*(\d+);", "bitsPerCell")
    v_rows = read_back(PARAM_CPP, r"numRowSubArray\s*=\s*(\d+);", "numRowSubArray")
    v_cols = read_back(PARAM_CPP, r"numColSubArray\s*=\s*(\d+);", "numColSubArray")
    print(f"   verified in Param.cpp: bitsPerCell={v_bits}, "
          f"numRowSubArray={v_rows}, numColSubArray={v_cols}")

    return logical_rows


def build():
    """make clean && make, with the return code actually checked."""
    subprocess.run("make clean", shell=True, cwd="NeuroSIM",
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    r = subprocess.run("make", shell=True, cwd="NeuroSIM",
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        tail = "\n".join(r.stdout.strip().splitlines()[-30:])
        raise RuntimeError(
            "make FAILED -- refusing to run.\n"
            "If this is ignored, the old binary stays in place and the run\n"
            "reports the PREVIOUS config's numbers, which look perfectly "
            "plausible.\n--- last 30 lines of make output ---\n" + tail)
    print("   build ok")


def run_inference(tag):
    """Run inference.py and interrupt it once the first HW batch is done."""
    cmd = INFER_CMD.format(model=MODEL_NAME)
    infer_log = os.path.join(LOG_DIR, f"inference_{tag}.log")

    env = dict(os.environ); env["BLOCKSTATS_TAG"] = tag
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            preexec_fn=os.setsid, env=env)

    hw_started = False
    hit_mark = False
    t0 = time.time()

    try:
        with open(infer_log, "w") as lf:
            for line in proc.stdout:
                lf.write(line)
                lf.flush()
                sys.stdout.write(line)

                if "Running inference with detailed hardware simulation" in line:
                    hw_started = True

                if hw_started and ("4%|" in line or "1/25" in line
                                   or "800 test images" in line):
                    hit_mark = True
                    print("\n-> [auto-interrupt] first HW batch done; "
                          "letting CSV writes flush...")
                    time.sleep(FLUSH_WAIT_S)
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                    break

                if time.time() - t0 > INFER_TIMEOUT_S:
                    print("\n-> [watchdog] inference exceeded the timeout; killing it")
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                    break
    finally:
        try:
            proc.wait(timeout=60)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            proc.wait()

    if not hit_mark:
        print("   WARNING: never saw the batch-completion marker. The layer "
              "CSVs may be incomplete.")
    return hit_mark


def run_trace(tag):
    """Run the NeuroSim hardware trace and save the log."""
    trace_sh = f"./layer_record_{MODEL_NAME}/trace_command.sh"
    if not os.path.exists(trace_sh):
        raise RuntimeError(f"{trace_sh} not found -- inference.py did not "
                           f"produce the layer records for this config")

    log_path = os.path.join(LOG_DIR, f"results_{MODEL_NAME}_{tag}.txt")
    env = dict(os.environ); env["BLOCKSTATS_TAG"] = tag   # labels the block report
    with open(log_path, "w") as lf:
        subprocess.run(f"bash {trace_sh}", shell=True,
                       stdout=lf, stderr=subprocess.STDOUT, env=env)

    size = os.path.getsize(log_path)
    with open(log_path) as f:
        head = f.read(4000)
    for marker, why in (
            ("SubArray Size is too large", "chip hierarchy guard fired"),
            ("Segmentation fault",          "main segfaulted"),
            ("core dumped",                 "main dumped core")):
        if marker in head:
            raise RuntimeError(
                f"{why}: {log_path} contains \"{marker}\". "
                f"The run produced NO usable numbers.")
    if size < 512:
        print(f"   WARNING: {log_path} is only {size} B -- the trace probably failed")
    return log_path, size


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="apply and verify the source edits only; no build, no run")
    ap.add_argument("--only", action="append", default=[],
                    help="restrict to a geometry, e.g. --only 512x512 (repeatable)")
    ap.add_argument("--n", type=int, action="append", default=[],
                    help="restrict to an n value (repeatable)")
    ap.add_argument("--no-restore", action="store_true",
                    help="leave the sources on the last config instead of restoring")
    args = ap.parse_args()

    geoms = GEOMETRIES
    if args.only:
        want = {s.lower().replace(" ", "") for s in args.only}
        geoms = [g for g in GEOMETRIES if f"{g[0]}x{g[1]}" in want]
        if not geoms:
            sys.exit(f"--only matched nothing. Available: "
                     f"{', '.join(f'{r}x{c}' for r, c in GEOMETRIES)}")
    ns = args.n or N_VALUES

    os.makedirs(LOG_DIR, exist_ok=True)
    configs = [(n, r, c) for n in ns for (r, c) in geoms]

    print(f"\n{len(configs)} config(s): "
          + ", ".join(f"n={n} {r}x{c}" for n, r, c in configs))
    print(f"rows are treated as {'PHYSICAL' if ROWS_ARE_PHYSICAL else 'LOGICAL'}\n")

    saved = backup_sources()
    results = []

    try:
        for i, (n, rows, cols) in enumerate(configs, 1):
            tag = f"{CELL_LABEL.format(n=n)}_{rows}x{cols}"
            log_path = os.path.join(LOG_DIR, f"results_{MODEL_NAME}_{tag}.txt")

            print("=" * 72)
            print(f"[{i}/{len(configs)}]  n={n}  {rows} x {cols}   (tag {tag})")
            print("=" * 72)

            if os.path.exists(log_path) and not args.dry_run:
                print(f"-> already done ({log_path}); skipping. "
                      f"Delete it to re-run.")
                results.append((tag, "skipped", os.path.getsize(log_path)))
                continue

            try:
                print("-> applying config")
                apply_config(n, rows, cols)

                if args.dry_run:
                    results.append((tag, "edits ok (dry run)", 0))
                    continue

                print("-> building")
                build()

                print("-> running inference")
                ok = run_inference(tag)

                print("-> running hardware trace")
                path, size = run_trace(tag)
                results.append((tag, "ok" if ok else "ok (marker missed)", size))
                print(f"-> saved {path} ({size} B)\n")

            except (EditError, RuntimeError) as e:
                print(f"\n!! config {tag} FAILED:\n{e}\n")
                results.append((tag, f"FAILED: {str(e).splitlines()[0]}", 0))
                # An edit failure means every later config is suspect too.
                if isinstance(e, EditError):
                    print("!! aborting the sweep -- a source edit did not apply, "
                          "so later runs would silently reuse stale values.")
                    break

    except KeyboardInterrupt:
        print("\n-> interrupted by user")
    finally:
        if not args.no_restore:
            restore_sources(saved)
        else:
            print("-> --no-restore: sources left on the last applied config")

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    print(f"{'config':24s} {'status':34s} {'log bytes':>10s}")
    for tag, status, size in results:
        print(f"{tag:24s} {status:34s} {size:10d}")
    print(f"\nLogs in ./{LOG_DIR}/")


if __name__ == "__main__":
    main()
