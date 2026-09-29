#!/usr/bin/env bash
# sweep_paper.sh -- one-dimensional sweeps for paper figures.
#
# Each FAMILY varies exactly ONE thing and holds everything else fixed, so every
# curve has an unambiguous x axis.  The grid sweep in sweep.sh is for auditing;
# this one is for plotting.
#
#   square   rows = cols, n fixed            -> "how does a square subarray scale?"
#   rows     rows varies, cols and n fixed   -> isolates the COLUMN-direction lines
#   cols     cols varies, rows and n fixed   -> isolates the ROW-direction lines
#   planes   n varies, rows = cols fixed     -> the 3D stacking result
#   aspect   rows x cols held CONSTANT       -> shape at fixed capacity
#
# Write scheme: VANILLA only (scheme 3) -- the standard 1T1C single plate cycle.
#
# Usage:
#   ./sweep_paper.sh                  # 2T-nC qndro, all five families
#   ./sweep_paper.sh 5 2 planes       # just one family
#   ./sweep_paper.sh 6 3              # 1T-nC dro, all five
#
# Outputs (one CSV per family, same 143-column schema as sweep.sh):
#   paper_<cell>_<rdo>_square.csv   ... _rows.csv  _cols.csv  _planes.csv  _aspect.csv
#   paper_<cell>_<rdo>_logs/        full text report per run

set -u
CELL=${1:-5}
RDO=${2:-2}
ONLY=${3:-all}                 # all | square | rows | cols | planes | aspect

VMIN=0.025                     # offset-cancelled sense amp; 0 = use Param's value
SCHEME=3                       # 3 = vanilla.  2 = erase-program, 1 = differential
P1=0.5                         # fraction of 1s already stored
Q1=0.5                         # fraction of 1s being written

PREFIX=paper_${CELL}_${RDO}
LOGS=${PREFIX}_logs
mkdir -p "$LOGS"

[ -x ./rw_main ] || { echo "build it first:  make rw_main"; exit 1; }

# --------------------------------------------------------------- sweep axes --
# Start at 16, not 8: below ~16 rows NeuroSim's area model returns large negative
# decoder areas (the array is narrower than one peripheral cell row).  16 and 32
# are included so you can SEE where it breaks -- plot_paper.py drops any point
# whose macro area is non-positive from the area and density figures and tells
# you which ones it dropped.
SQUARE_SIDES="16 32 64 128 256 512 1024"
SQUARE_N=64

ROWS_LIST="16 32 64 128 256 512 1024"
ROWS_COLS=256                  # held fixed while rows varies
ROWS_N=64

COLS_LIST="16 32 64 128 256 512 1024"
COLS_ROWS=256                  # held fixed while cols varies
COLS_N=64

PLANES_LIST="2 4 8 16 32 64 128 256 512"
PLANES_ROWS=256
PLANES_COLS=256

# constant capacity 256 x 256 = 65536 cells per plane, reshaped
ASPECT_PAIRS="16:4096 32:2048 64:1024 128:512 256:256 512:128 1024:64 2048:32 4096:16"
ASPECT_N=64

BLOCKS="rblDecoder rblSwitchMatrix rslDecoder rslSwitchMatrix \
wblDecoder wblSwitchMatrix wwlDecoder wwlSwitchMatrix \
wplDecoder wplSwitchMatrix sslDecoder sslSwitchMatrix \
currentSenseAmp outputLatch cellArray"

# ---------------------------------------------------------------- header -----
# IMPORTANT: must stay byte-identical to sweep.sh's header, or build_rw_xlsx.py
# and plot_paper.py will mis-parse.  If you add a block, add it in BOTH scripts.
emit_header() {
  {
    printf 'cell,rdo,n,rows,cols,p1,q1,scheme,flipFrac'
    printf ',eReadCore_J,eReadMacro_J,eReadMacro_fJb'
    printf ',eWriteCore_J,eWriteMacro_J,eWriteMacro_fJb'
    printf ',tReadCore_s,tReadMacro_s,tWriteCore_s,tWriteMacro_s'
    printf ',areaCore_um2,areaMacro_um2,leakMacro_W'
    printf ',dVsense_V,senseVerdict'
    printf ',eDecRd_J,eDecWr_J,eSenseAmp_J,ePrecharge_J,eLatch_J'
    printf ',capFE_F,capNode_F,Qsw_C,vIntBias_V,vUnsel_V,iRead1_A,iRead0_A'
    printf ',tauSwitch_s,numWritePulses'
    printf ',eRd_cellSwitch_J,eRd_cellLinear_J,eRd_transistorDC_J,eRd_lines_J,eRd_restore_J'
    printf ',eWr_cells_J,eWr_inhibit_J,eWr_lines_J'
    printf ',driveFrac,switchesPerWrite,writesToFailure,eWriteTrue_J,tWriteTrue_s'
    printf ',t_select_s,t_charge_s,t_settle_s,t_senseAmp_s,t_latch_s,t_restore_s,t_selectWrite_s'
    for b in $BLOCKS; do
        printf ',%s_eRead_J,%s_eWrite_J,%s_tRead_s,%s_tWrite_s,%s_leak_W,%s_area_m2' \
               "$b" "$b" "$b" "$b" "$b" "$b"
    done
    printf '\n'
  } > "$1"
}

# ----------------------------------------------------------------- one run ---
NFAIL=0
run_one() {   # $1=n  $2=rows  $3=cols  $4=csv  $5=family
  local N=$1 R=$2 C=$3 CSV=$4 FAM=$5
  local TAG="${FAM}_n${N}_${R}x${C}"
  ./rw_main "$CELL" "$RDO" "$N" "$R" "$C" "$P1" "$Q1" "$SCHEME" 0 "$VMIN" \
       > "${LOGS}/${TAG}.txt" 2>> "$CSV"
  local rc=$?
  if [ $rc -gt 2 ]; then
      echo "    FAILED rc=$rc  $TAG"
      NFAIL=$((NFAIL+1))
  elif [ $rc -eq 2 ]; then
      echo "    sense FAIL (kept, plotted as an open marker)  $TAG"
  fi
  printf '.'
}

want() { [ "$ONLY" = "all" ] || [ "$ONLY" = "$1" ]; }

# --------------------------------------------------------------- families ----
if want square; then
  CSV=${PREFIX}_square.csv; emit_header "$CSV"
  echo "square   n=${SQUARE_N}, rows=cols in { $SQUARE_SIDES }"
  for s in $SQUARE_SIDES; do run_one "$SQUARE_N" "$s" "$s" "$CSV" square; done
  echo " -> $CSV"
fi

if want rows; then
  CSV=${PREFIX}_rows.csv; emit_header "$CSV"
  echo "rows     n=${ROWS_N}, cols=${ROWS_COLS}, rows in { $ROWS_LIST }"
  for r in $ROWS_LIST; do run_one "$ROWS_N" "$r" "$ROWS_COLS" "$CSV" rows; done
  echo " -> $CSV"
fi

if want cols; then
  CSV=${PREFIX}_cols.csv; emit_header "$CSV"
  echo "cols     n=${COLS_N}, rows=${COLS_ROWS}, cols in { $COLS_LIST }"
  for c in $COLS_LIST; do run_one "$COLS_N" "$COLS_ROWS" "$c" "$CSV" cols; done
  echo " -> $CSV"
fi

if want planes; then
  CSV=${PREFIX}_planes.csv; emit_header "$CSV"
  echo "planes   rows=cols=${PLANES_ROWS}, n in { $PLANES_LIST }"
  for n in $PLANES_LIST; do run_one "$n" "$PLANES_ROWS" "$PLANES_COLS" "$CSV" planes; done
  echo " -> $CSV"
fi

if want aspect; then
  CSV=${PREFIX}_aspect.csv; emit_header "$CSV"
  echo "aspect   n=${ASPECT_N}, rows x cols = 65536 held constant"
  for p in $ASPECT_PAIRS; do
    run_one "$ASPECT_N" "${p%%:*}" "${p##*:}" "$CSV" aspect
  done
  echo " -> $CSV"
fi

echo
echo "done.  $NFAIL hard failures.  Full reports in ${LOGS}/"
echo
echo "next:"
echo "    python3 plot_paper.py ${PREFIX} -o figs"
echo "    python3 build_rw_xlsx.py ${PREFIX}_planes.csv -o ${PREFIX}_planes.xlsx"
