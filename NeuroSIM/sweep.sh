#!/usr/bin/env bash
# sweep.sh -- read/write operator mode sweep, with per-block detail.
#
#   n (stacked planes)       : 64, 256      (+ any explicit configs in EXTRA)
#   physical rows x columns  : 64x64, 64x512, 512x64, 512x512
#   write schemes            : 2 = erase-then-program   (blind, 2 phases)
#                              3 = vanilla single cycle (blind, 2 edges)  <-- baseline
#                              1 = differential RMW     (owes a read)
#
# Single-row read and single-row write throughout: rw_main sets
# activityRowRead = activityRowWrite = 1/numRow, so numRow x activity == 1.
#
# Usage:
#   ./sweep.sh                 # 2T-nC qndro, random data, all three write schemes
#   ./sweep.sh 6 3             # 1T-nC dro
#   ./sweep.sh 5 2 corners     # + the all-0 / all-1 data corners
#   ./sweep.sh 5 2 rows        # column-path scan: rows 64 -> 8192 at 512 columns
#                              #   (output rw_5_2_rows.csv, vanilla only)

set -u
CELL=${1:-5}
RDO=${2:-2}
MODE=${3:-plain}
VMIN=0.025                     # offset-cancelled sense amp; 0 = use Param's value
OUT=rw_${CELL}_${RDO}
[ "$MODE" = "rows" ] && OUT=${OUT}_rows
LOGS=${OUT}_logs
mkdir -p "$LOGS"

[ -x ./rw_main ] || { echo "build it first:  make rw_main"; exit 1; }

NS="64 256"
GEOMS="64:64 64:512 512:64 512:512"
EXTRA="8:8:8"                  # explicit n:rows:cols runs, appended to the grid
                               # e.g. EXTRA="8:8:8 16:128:128 512:256:256"
if [ "$MODE" = "rows" ]; then
    # SENSE/LINE FIX (2026-09-25): the row-count scan.  The column line is the
    # read transistor's source line, one amplifier per column: its capacitance
    # (t_col, t_dev) and its resistance (far-row IR degeneration) both scale
    # with numRow.  512 columns so the plate strip and the RBL page current are
    # the same in every run.
    GEOMS="64:512 128:512 256:512 512:512 1024:512 2048:512 4096:512 8192:512"
    EXTRA=""
fi

# FIX 4: planeDecoder MUST sit immediately after sslSwitchMatrix -- it is the
# last entry of rw_main's P[] array, and the emit loop runs over P[] before
# appending the sense amp / latch / cellArray groups by hand.  The same order
# is hardcoded in build_rw_xlsx.py's BLOCKS list.
BLOCKS="rblDecoder rblSwitchMatrix rslDecoder rslSwitchMatrix \
wblDecoder wblSwitchMatrix wwlDecoder wwlSwitchMatrix \
wplDecoder wplSwitchMatrix sslDecoder sslSwitchMatrix \
planeDecoder \
currentSenseAmp outputLatch cellArray"

if [ "$MODE" = "corners" ]; then
    PATTERNS="0.0:0.0 0.0:1.0 0.5:0.5 1.0:0.0 1.0:1.0"
else
    PATTERNS="0.5:0.5"
fi
SCHEMES="2 3 1"                # erase-program, vanilla, differential
[ "$MODE" = "rows" ] && SCHEMES="3"

# ---------------------------------------------------------------- header ----
{
  printf 'cell,rdo,n,rows,cols,p1,q1,scheme,flipFrac'
  printf ',eReadCore_J,eReadMacro_J,eReadMacro_fJb'
  printf ',eWriteCore_J,eWriteMacro_J,eWriteMacro_fJb'
  printf ',tReadCore_s,tReadMacro_s,tWriteCore_s,tWriteMacro_s'
  printf ',areaCore_um2,areaMacro_um2,leakMacro_W'
  printf ',areaArray_um2,areaStaircase_um2'
  printf ',dVsense_V,senseVerdict'
  printf ',eDecRd_J,eDecWr_J,eSenseAmp_J,ePrecharge_J,eLatch_J'
  # physics
  printf ',capFE_F,capNode_F,Qsw_C,vIntBias_V,vUnsel_V,iRead1_A,iRead0_A'
  printf ',tauSwitch_s,numWritePulses'
  # read breakdown
  printf ',eRd_cellSwitch_J,eRd_cellLinear_J,eRd_transistorDC_J,eRd_lines_J,eRd_restore_J'
  # write breakdown + the columns that separate the three schemes
  #
  # FIX 1 -- ORDER MUST MATCH rw_main.cpp's emit EXACTLY.  nHalfSel is emitted
  # LAST in this group, not fourth.  Naming it fourth rotated every field after
  # it by one position (csv.DictReader zips names to positions), which is how
  # "true latency incl. owed read" ended up reporting nHalfSel in seconds and
  # "reversals per cell per write" reported the endurance constant 2e12.
  printf ',eWr_cells_J,eWr_inhibit_J,eWr_lines_J'
  printf ',driveFrac,switchesPerWrite,writesToFailure,eWriteTrue_J,tWriteTrue_s'
  printf ',nHalfSel'
  # FIX 7 -- disturb margin.  V_half/V_coercive > 1 means the V/2 scheme is
  # writing the half-selected cells instead of inhibiting them.
  printf ',vHalfSel_V,eHalfSel_MVcm,disturbRatio'
  # latency breakdown
  printf ',t_select_s,t_charge_s,t_settle_s,t_senseAmp_s,t_latch_s,t_restore_s,t_selectWrite_s'
  printf ',t_switch_s,t_plateRC_s,t_transfer_s,t_settleRow_s,t_settleCol_s'
  printf ',t_senseInt_s,iDiff_A,iLeakRows_A,capSenseNode_F,readOv_V'
  # per block
  for b in $BLOCKS; do
      printf ',%s_eRead_J,%s_eWrite_J,%s_tRead_s,%s_tWrite_s,%s_leak_W,%s_area_m2' \
             "$b" "$b" "$b" "$b" "$b" "$b"
  done
  # SENSE/LINE FIX (2026-09-25) -- diagnostics APPENDED after the block group so
  # every awk field number above is unchanged.  ORDER MUST MATCH rw_main's
  # trailing emit.  t_transfer_s (above) is now always 0: the R_plate*n*C_FE
  # term was removed from the read chain; t_senseAmp_s = develop + regenerate.
  printf ',t_senseDev_s,t_senseRegen_s,tauLatch_s,t_senseNeuroSim_s,t_plateSupply_s'
  printf ',rPlateDrv_ohm,rPlateEff_ohm,capJunctionTr_F,capWBLpar_F,capWBLwire_F'
  printf ',eRd_unselStrips_J,eWr_unselStrips_J,senseRefFrac'
  # column-path / sense-mode diagnostics (second appended group, same rule)
  printf ',senseMode,gClamp_S,capColSense_F,iRead1Ideal_A,iRead0Ideal_A,rRslPath_ohm'
  printf ',vSourceIR_V,degenLoss_frac,vDropRBL_V,rblDropRatio,kLineRc'
  printf '\n'
} > "${OUT}.csv"

# ------------------------------------------------------------- run list -----
# one flat list of n:rows:cols triples = the grid, plus the explicit extras
RUNS=""
for n in $NS; do
  for g in $GEOMS; do RUNS="$RUNS ${n}:${g}"; done
done
RUNS="$RUNS $EXTRA"

# ----------------------------------------------------------------- runs -----
n_runs=0
for r in $RUNS; do
  N=${r%%:*}; tmp=${r#*:}; ROWS=${tmp%%:*}; COLS=${tmp##*:}
  for s in $SCHEMES; do
    for p in $PATTERNS; do
      P1=${p%%:*}; Q1=${p##*:}
      TAG="n${N}_${ROWS}x${COLS}_p${P1}_q${Q1}_s${s}"
      ./rw_main "$CELL" "$RDO" "$N" "$ROWS" "$COLS" "$P1" "$Q1" "$s" 0 "$VMIN" \
           > "${LOGS}/${TAG}.txt" 2>> "${OUT}.csv"
      rc=$?
      [ $rc -gt 2 ] && echo "  FAILED rc=$rc  $TAG"
      n_runs=$((n_runs+1))
    done
  done
done

NCOL=$(head -1 "${OUT}.csv" | awk -F, '{print NF}')
echo
echo "$n_runs runs, $NCOL columns -> ${OUT}.csv   (full reports in ${LOGS}/)"
echo

# --------------------------------------------------------- headline table ---
awk -F, 'NR>1 && $6=="0.5" && $8=="erase-program" {
    printf "%5s %6s %6s  %10.3f %12.3f  %9.3f %10.3f  %11.1f %8.2f %6s\n",
           $3,$4,$5, $12, $15, $17*1e9, $19*1e9, $21, $25*1e3, $26
}' "${OUT}.csv" | sort -n | \
awk 'BEGIN{printf "%5s %6s %6s  %10s %12s  %9s %10s  %11s %8s %6s\n",
     "n","rows","cols","read fJ/b","write fJ/b","tRead ns","tWrite ns","macro um2","dV mV","sense";
     printf "  (write column is erase-then-program; see the scheme table below)\n"}
     {print}'

# ------------------------------------------------- write-scheme comparison ---
echo
printf '%5s %6s %6s | %23s | %17s | %s\n' \
       "n" "rows" "cols" "     write fJ/bit      " " switches/write  " " writes to failure"
printf '%5s %6s %6s | %7s %7s %7s | %5s %5s %5s | %9s %9s %9s\n' \
       "" "" "" "erase" "vanil" "diff" "erase" "vanil" "diff" "erase" "vanil" "diff"
awk -F, 'NR>1 && $6=="0.5" {
    k = sprintf("%05d %05d %05d", $3, $4, $5)
    wfb[k "|" $8] = $15          # eWriteMacro_fJb
    sw [k "|" $8] = $50          # switchesPerWrite   (rw_main emit field 50)
    wtf[k "|" $8] = $51          # writesToFailure    (rw_main emit field 51)
    seen[k] = 1
}
END{
    for (k in seen) {
        split(k, f, " ")
        printf "%5d %6d %6d | %7.2f %7.2f %7.2f | %5.3f %5.3f %5.3f | %9.2e %9.2e %9.2e\n",
               f[1], f[2], f[3],
               wfb[k"|erase-program"], wfb[k"|vanilla"], wfb[k"|differential"],
               sw [k"|erase-program"], sw [k"|vanilla"], sw [k"|differential"],
               wtf[k"|erase-program"], wtf[k"|vanilla"], wtf[k"|differential"]
    }
}' "${OUT}.csv" | sort -k1,1n -k2,2n -k3,3n

cat <<'NOTE'

  Energy barely separates the three -- the V/2 half-select dominates the write and
  all three take 2 plate phases. ENDURANCE is where they differ: erase-then-program
  switches 1.0 per write however little the data changes, while vanilla gets the
  Hamming count for free (no read) and differential pays a read for the same count.
NOTE

# ------------------------------------------------- disturb margin (FIX 7) ---
# V_half/V_coercive > 1 means the "inhibited" cells are above the coercive field,
# i.e. the V/2 scheme is writing them.  That is a design failure, not a rounding
# error, so it gets its own line rather than hiding in the inhibit energy.
echo
NF_ALL=$(head -1 "${OUT}.csv" | awk -F, '{print NF}')
IDX_HS=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^nHalfSel$'     | cut -d: -f1)
IDX_VH=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^vHalfSel_V$'   | cut -d: -f1)
IDX_DR=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^disturbRatio$' | cut -d: -f1)
awk -F, -v hs="$IDX_HS" -v vh="$IDX_VH" -v dr="$IDX_DR" '
NR==1 { printf "%5s %6s %6s | %11s %9s %9s  %s\n",
        "n","rows","cols","half-sel","V_half","V_h/V_c","verdict"; next }
$6=="0.5" && $8=="vanilla" {
    printf "%5d %6d %6d | %11d %9.3f %9.3f  %s\n",
           $3,$4,$5, $hs, $vh, $dr, ($dr > 1.0 ? "OVER COERCIVE -- not inhibiting" : "ok")
}' "${OUT}.csv" | sort -k1,1n -k2,2n -k3,3n -s

# ------------------------------------------------- column path (2026-09-25) ---
# The row count enters the read through the COLUMN line (the read transistor's
# source line): its capacitance sets the column settle + development time and
# its resistance degenerates the far row.  The page current through the ROW
# line (RBL) is a column-count effect; a ratio > 1 means the far columns lose
# V_DS and a page-wide current read is not viable at this RBL width.
echo
IDX_TC=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^t_settleCol_s$'   | cut -d: -f1)
IDX_TD=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^t_senseDev_s$'    | cut -d: -f1)
IDX_TA=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^t_senseAmp_s$'    | cut -d: -f1)
IDX_CC=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^capColSense_F$'   | cut -d: -f1)
IDX_DL=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^degenLoss_frac$'  | cut -d: -f1)
IDX_VR=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^rblDropRatio$'    | cut -d: -f1)
IDX_TR=$(head -1 "${OUT}.csv" | tr ',' '\n' | grep -n '^tReadMacro_s$'    | cut -d: -f1)
awk -F, -v tc="$IDX_TC" -v td="$IDX_TD" -v ta="$IDX_TA" -v cc="$IDX_CC" -v dl="$IDX_DL" -v vr="$IDX_VR" -v tr="$IDX_TR" '
NR==1 { printf "%5s %6s %6s | %8s %10s %10s %10s | %9s %9s | %10s  %s\n",
        "n","rows","cols","C_col fF","t_col ns","t_dev ns","t_sense ns","far-row","RBL drop","tRead ns","verdict"; next }
$6=="0.5" && $8=="vanilla" {
    printf "%5d %6d %6d | %8.2f %10.4f %10.4f %10.4f | %8.2f%% %9.3f | %10.3f  %s\n",
           $3,$4,$5, $cc*1e15, $tc*1e9, $td*1e9, $ta*1e9, 100*$dl, $vr, $tr*1e9,
           ($vr > 1.0 ? "RBL IR: page current collapses V_DS" : "ok")
}' "${OUT}.csv" | sort -k1,1n -k2,2n -k3,3n -s
