#!/usr/bin/env python3
"""
apply_blockstats.py -- add chip-level per-block latency/energy accounting

WHAT IT DOES
    1. copies BlockStats.h into NeuroSIM/
    2. adds  #include "BlockStats.h"  to SubArray.cpp
    3. inserts ONE accumulator block at the end of the 2T-nC branch of
       SubArray::CalculatePower

    Nothing is added to main.cpp: the accumulator's destructor prints the
    report at program exit.

WHY ONE HOOK IS ENOUGH
    ProcessingUnit.cpp:370-378 calls
        subArray->CalculateLatency(...);      <- fills every .readLatency
        if (!CalculateclkFreq) {
            subArray->CalculatePower(...);    <- fills every .readDynamicEnergy
            *readDynamicEnergy += subArray->readDynamicEnergy;
        }
    on the same object, inside the loop over input vectors, which sits inside
    the loops over subarrays (ProcessingUnit), PEs (Tile.cpp:479) and tiles
    (Chip.cpp:765).  Accumulating at the end of CalculatePower therefore sees
    exactly the same set of contributions the chip total is built from, and
    only on pass 2.

READING THE REPORT
    ENERGY per block sums to the chip read/write energy -- it is additive all
    the way up.
    LATENCY per block is total block-time over every invocation.  Tile and Chip
    combine latencies with MAX(), not +, so this is NOT the critical path; the
    report prints the MAX separately so the gap is visible.
    LEAKAGE is a power, tracked as a MAX (per-subarray value), not a sum.

    Set BLOCKSTATS_TAG in the environment to label each run, e.g.
        BLOCKSTATS_TAG="n=64 512x512" ./NeuroSIM/main ...

USAGE
    python3 apply_blockstats.py /path/to/NeuroSIM/SubArray.cpp
    (BlockStats.h must sit next to this script)
"""
import os, shutil, sys

HOOK = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "hook.inc")).read().rstrip("\n")

ANCHOR = "                \treadDynamicEnergyOther = readDynamicEnergy - readDynamicEnergyADC - readDynamicEnergyAccum;"

INC_OLD = '#include "SubArray.h"'
INC_NEW = '#include "SubArray.h"\n#include "BlockStats.h"'


def main():
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(1)
    sa = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    hdr_src = os.path.join(here, "BlockStats.h")
    if not os.path.exists(hdr_src):
        sys.exit("BlockStats.h not found next to this script")
    if not os.path.exists(sa):
        sys.exit(f"{sa} not found")

    raw = open(sa, newline="").read()
    crlf = "\r\n" in raw
    src = raw.replace("\r\n", "\n") if crlf else raw

    if "BlockStats.h" in src:
        sys.exit("SubArray.cpp already includes BlockStats.h -- already patched. "
                 "Restore from SubArray.cpp.bak_blockstats first.")

    n_anchor = src.count(ANCHOR)
    if n_anchor == 0:
        sys.exit("FAIL: could not find the end-of-2T-nC-power anchor:\n  " + ANCHOR)
    if src.count(INC_OLD) != 1:
        sys.exit(f"FAIL: expected exactly one '{INC_OLD}', found {src.count(INC_OLD)}")

    # header
    dest = os.path.join(os.path.dirname(os.path.abspath(sa)), "BlockStats.h")
    shutil.copyfile(hdr_src, dest)
    print(f"  ok    copied BlockStats.h -> {dest}")

    shutil.copyfile(sa, sa + ".bak_blockstats")

    src = src.replace(INC_OLD, INC_NEW, 1)
    print("  ok    added #include \"BlockStats.h\"")

    src = src.replace(ANCHOR, ANCHOR + "\n\n" + HOOK, 1)
    print(f"  ok    inserted accumulator at the FIRST of {n_anchor} anchor(s) "
          f"(the 2T-nC branch)")

    if crlf:
        src = src.replace("\n", "\r\n")
    open(sa, "w", newline="").write(src)
    print(f"\nDone. Backup: {os.path.basename(sa)}.bak_blockstats")
    print("Rebuild with:  cd NeuroSIM && make clean && make")
    print("The report prints itself at the end of every ./NeuroSIM/main run.")


if __name__ == "__main__":
    main()
