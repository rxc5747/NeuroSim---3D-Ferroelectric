#!/usr/bin/env bash
# sweep_all.sh -- 3D FeRAM, SRAM and RRAM into ONE csv, in ONE set of tiers,
#                 for the technology-comparison figures.
#
# Everything goes through rw_main so the tier definitions are identical:
#   CORE  = cells + the array lines they charge
#   MACRO = core + the SAME readout periphery for every technology
#           (one current sense amp per column + one output latch)
#           + that technology's own address path
#
# The baselines are 2D, so n = 1 and there is no write-scheme structure --
# one scheme each.  The FeRAM runs keep their plane and scheme sweeps.
set -u
OUT=compare
LOGS=${OUT}_logs
VMIN=0.025
mkdir -p "$LOGS"
[ -x ./rw_main ] || { echo "build it first:  make clean && make rw_main"; exit 1; }

GEOMS="64:64 64:512 512:64 512:512"

#  label          cell rdo  n-list      schemes
CELLS="
2T-nC-qndro         5   2   64,256      2,3,1
1T-nC-dro           6   3   64,256      3
SRAM-6T             1   3   1           3
RRAM-1T1R           2   3   1           3
"

# ---- header: reuse sweep.sh's, so build_rw_xlsx.py is unchanged ------------
if [ -f rw_5_2.csv ]; then
    head -1 rw_5_2.csv > "${OUT}.csv"
else
    echo "run  ./sweep.sh 5 2  once first so the CSV header exists"; exit 1
fi
NHDR=$(head -1 "${OUT}.csv" | awk -F, '{print NF}')

n_runs=0
while read -r LABEL CELL RDO NLIST SLIST; do
  [ -z "${LABEL:-}" ] && continue
  for N in $(echo "$NLIST" | tr ',' ' '); do
    for g in $GEOMS; do
      ROWS=${g%%:*}; COLS=${g##*:}
      for s in $(echo "$SLIST" | tr ',' ' '); do
        TAG="${LABEL}_n${N}_${ROWS}x${COLS}_s${s}"
        ./rw_main "$CELL" "$RDO" "$N" "$ROWS" "$COLS" 0.5 0.5 "$s" 0 "$VMIN" \
             > "${LOGS}/${TAG}.txt" 2>> "${OUT}.csv"
        rc=$?
        [ $rc -gt 2 ] && echo "  FAILED rc=$rc  $TAG"
        n_runs=$((n_runs+1))
      done
    done
  done
done <<< "$CELLS"

NDAT=$(sed -n 2p "${OUT}.csv" | awk -F, '{print NF}')
echo
echo "$n_runs runs -> ${OUT}.csv"
echo "header fields: $NHDR      data fields: $NDAT"
[ "$NHDR" != "$NDAT" ] && echo "  *** MISMATCH -- rebuild rw_main, the binary is stale ***"

# ---- the guard that matters ------------------------------------------------
echo
echo "zero-tier check (any row here is a commented-out path, NOT a result):"
awk -F, 'NR>1 && (($10+0)==0 || ($13+0)==0 || ($17+0)==0 || ($19+0)==0) {
    printf "  %-26s n=%-4s %4sx%-4s  eRdCore=%s eWrCore=%s tRd=%s tWr=%s\n",$1,$3,$4,$5,$10,$13,$17,$19 }' \
    "${OUT}.csv" | grep . || echo "  none -- every technology produced a real number"

# ---- headline comparison ---------------------------------------------------
echo
printf '%-26s %5s %6s %6s | %10s %11s | %9s %10s | %12s\n' \
       "technology" "n" "rows" "cols" "read fJ/b" "write fJ/b" "tRead ns" "tWrite ns" "core um^2"
awk -F, 'NR>1 && $8=="vanilla" {
    printf "%-26s %5s %6s %6s | %10.4f %11.3f | %9.4f %10.4f | %12.2f\n",
           $1,$3,$4,$5,$12,$15,$17*1e9,$19*1e9,$20 }' "${OUT}.csv" | sort -t'|' -k1,1
