#!/usr/bin/env python3
"""
rw_sweep.py -- read/write operator mode sweep over the 8-config grid.

    n (stacked planes)      : 64, 256
    physical rows x columns : 64x64, 64x512, 512x64, 512x512

Unlike the CIM sweep this needs NO netStructure file, NO Chip/Tile/PE hierarchy
and NO Param.cpp editing, so none of the "SubArray Size is too large" /segfault
failure modes apply.  Every run is one process, a few milliseconds.

For each geometry it runs the data-pattern corners as well as the random point,
because the whole point of this mode is that the energy is data dependent:

    read  : p1 = 0.0 (all zeros), 0.5 (random), 1.0 (all ones)
    write : (p1,q1) = (0,0) (0,1) (1,0) (1,1) (0.5,0.5)
    scheme: 2 = erase-then-program, 1 = differential

Usage
    ./rw_sweep.py                       # full sweep, default binary ./rw_main
    ./rw_sweep.py --bin ../NeuroSIM/rw_main --out rw
    ./rw_sweep.py --cell 6 --rdo 3      # 1T-nC DRO instead of 2T-nC QNDRO
"""

import argparse, csv, itertools, os, subprocess, sys

CSV_HEADER = [
    "cell", "rdo", "n", "rows", "cols", "p1", "q1", "scheme", "flipFrac",
    "eReadCore_J", "eReadMacro_J", "eReadMacro_fJ_per_bit",
    "eWriteCore_J", "eWriteMacro_J", "eWriteMacro_fJ_per_bit",
    "tReadCore_s", "tReadMacro_s", "tWriteCore_s", "tWriteMacro_s",
    "areaCore_um2", "areaMacro_um2", "leakMacro_W",
    "dVsense_V", "senseVerdict",
    "eDecodeRead_J", "eDecodeWrite_J", "eSenseAmp_J", "ePrecharge_J", "eLatch_J",
]

GEOMS = [(n, r, c)
         for n in (64, 256)
         for (r, c) in ((64, 64), (64, 512), (512, 64), (512, 512))]

# (p1, q1) corners.  p1 also drives the read, so a single run gives both.
PATTERNS = [(0.0, 0.0), (0.0, 1.0), (0.5, 0.5), (1.0, 0.0), (1.0, 1.0)]


def run(binary, cell, rdo, n, rows, cols, p1, q1, scheme, logdir):
    tag = f"{cell}_{rdo}_n{n}_{rows}x{cols}_p{p1}_q{q1}_s{scheme}"
    cmd = [binary, str(cell), str(rdo), str(n), str(rows), str(cols),
           str(p1), str(q1), str(scheme)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except FileNotFoundError:
        sys.exit(f"binary not found: {binary}\n"
                 f"build it with:  make rw_main   (see RW_MODE_GUIDE.md part 6)")
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT  {tag}")
        return None
    if logdir:
        with open(os.path.join(logdir, tag + ".txt"), "w") as f:
            f.write(r.stdout)
    # rc 2 means the sense-margin DRC failed; the row is still emitted, flagged.
    if r.returncode not in (0, 2):
        print(f"  FAILED rc={r.returncode}  {tag}")
        if r.stderr.strip():
            print("   ", r.stderr.strip().splitlines()[-1])
        return None
    line = [l for l in r.stderr.strip().splitlines() if l.count(",") >= 20]
    if not line:
        print(f"  NO CSV  {tag}")
        return None
    return line[-1].split(",")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="./rw_main")
    ap.add_argument("--out", default="rw_sweep")
    ap.add_argument("--cell", type=int, default=5, help="5=2T-nC 6=1T-nC 7=1T1C")
    ap.add_argument("--rdo", type=int, default=0,
                    help="1=ndro 2=qndro 3=dro; 0 = pick the legal one for the cell")
    ap.add_argument("--schemes", default="2,1", help="write schemes to sweep")
    ap.add_argument("--logs", default="rw_logs", help="'' to skip full reports")
    args = ap.parse_args()

    # 1T-nC is DRO-only, 2T-nC is QNDRO.  Enforce it rather than letting a
    # meaningless combination into the csv.
    rdo = args.rdo or (3 if args.cell == 6 else 2)
    if args.cell == 6 and rdo != 3:
        sys.exit("1T-nC is DRO only (rdo=3)")
    if args.cell == 5 and rdo not in (2, 3):
        print("note: 2T-nC is normally QNDRO (rdo=2)")

    logdir = args.logs
    if logdir:
        os.makedirs(logdir, exist_ok=True)

    schemes = [int(s) for s in args.schemes.split(",")]
    rows_out = []
    print(f"sweeping cell={args.cell} rdo={rdo} over "
          f"{len(GEOMS)} geometries x {len(PATTERNS)} patterns x {len(schemes)} schemes "
          f"= {len(GEOMS)*len(PATTERNS)*len(schemes)} runs\n")

    for (n, r, c) in GEOMS:
        for scheme in schemes:
            for (p1, q1) in PATTERNS:
                out = run(args.bin, args.cell, rdo, n, r, c, p1, q1, scheme, logdir)
                if out:
                    rows_out.append(out[:len(CSV_HEADER)])

    csv_path = args.out + ".csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(CSV_HEADER)
        w.writerows(rows_out)
    print(f"\nwrote {csv_path}  ({len(rows_out)} rows)")

    # ---- headline table: the random-data point, erase-then-program ----------
    def f(x):
        try:
            return float(x)
        except Exception:
            return float("nan")

    print(f"\n{'n':>5}{'rows':>7}{'cols':>7}"
          f"{'read fJ/bit':>14}{'write fJ/bit':>14}"
          f"{'tRead ns':>11}{'tWrite ns':>11}"
          f"{'macro um2':>13}{'dV mV':>9}{'sense':>7}")
    for row in rows_out:
        if row[5] != "0.5" or row[7] != "erase-program":
            continue
        print(f"{row[2]:>5}{row[3]:>7}{row[4]:>7}"
              f"{f(row[11]):>14.3f}{f(row[14]):>14.3f}"
              f"{f(row[16])*1e9:>11.3f}{f(row[18])*1e9:>11.3f}"
              f"{f(row[20]):>13.1f}{f(row[22])*1e3:>9.2f}{row[23]:>7}")

    # ---- data dependence, at the largest geometry --------------------------
    print("\ndata dependence at the largest geometry (n=256, 512x512, "
          "erase-then-program):")
    print(f"  {'p1':>5}{'q1':>5}{'flip':>8}{'read fJ/bit':>14}{'write fJ/bit':>14}")
    for row in rows_out:
        if row[2] != "256" or row[3] != "512" or row[4] != "512":
            continue
        if row[7] != "erase-program":
            continue
        print(f"  {f(row[5]):>5.2f}{f(row[6]):>5.2f}{f(row[8]):>8.3f}"
              f"{f(row[11]):>14.3f}{f(row[14]):>14.3f}")


if __name__ == "__main__":
    main()
