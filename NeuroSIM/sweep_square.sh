#!/usr/bin/env bash
# sweep_square.sh -- read/write operator mode, SQUARE subarrays only.
#
# The two design axes of the paper figures:
#   subarray size  S : rows == cols == S        (default 32 64 128 256 512)
#   stacked planes n :                          (default 8 16 32 64 128 256 512)
# for each of the three write schemes
#   2 = erase-then-program (blind, 2 phases)
#   3 = vanilla single cycle (blind, 2 edges)   <-- baseline
#   1 = differential RMW   (owes a read)
#
# Single-row read and single-row write throughout: rw_main sets
# activityRowRead = activityRowWrite = 1/numRow, so numRow x activity == 1.
#
# The CSV header below is IDENTICAL to sweep.sh's (same rw_main emit order),
# so build_rw_xlsx.py and plot_rw_sweep.py read either file.
#
# Usage:
#   ./sweep_square.sh                    # 2T-nC qndro, random data, all schemes
#   ./sweep_square.sh 6 3                # 1T-nC dro
#   ./sweep_square.sh 5 2 corners        # + the all-0 / all-1 data corners
#
# Every axis is overridable from the environment without editing the file:
#   SIZES="64 128 256" NS="16 64 256" ./sweep_square.sh
#   SCHEMES="3"                         ./sweep_square.sh     # vanilla only
#   EXTRA="8:8:8 64:64:512"             ./sweep_square.sh     # extra n:rows:cols
#   VMIN=0.020                          ./sweep_square.sh
#
# Then:
#   python3 build_rw_xlsx.py rw_sq_5_2.csv
#   python3 plot_rw_sweep.py rw_sq_5_2.csv --outdir figs

set -u
CELL=${1:-5}
RDO=${2:-2}
MODE=${3:-plain}
VMIN=${VMIN:-0.025}            # offset-cancelled sense amp; 0 = use Param's value
OUT=${OUT:-rw_sq_${CELL}_${RDO}}
LOGS=${OUT}_logs
mkdir -p "$LOGS"

[ -x ./rw_main ] || { echo "build it first:  make rw_main"; exit 1; }

# ------------------------------------------------------------- the grid -----
# Sizes below 32 are kept out of the default grid on purpose: under ~16 rows
# NeuroSim's area model returns negative decoder areas (the 8x8x8 run in the
# old sweep), and 16 is borderline.  Add them through SIZES= or EXTRA= if you
# want the sense/energy numbers anyway -- build_rw_xlsx.py flags the area.
#
# n = 512 is expected to FAIL the sense check at VMIN = 25 mV on a flat pillar
# (dV falls as 1/n: 117 mV at n=64, 29 mV at n=256, ~15 mV at n=512).  It is in
# the default grid because "where the flat pillar stops being readable" is a
# result worth plotting, not a run to hide.  rw_main exits with rc <= 2 on a
# sense FAIL and the row is still written, so the figure gets the point.
SIZES=${SIZES:-"32 64 128 256 512"}
NS=${NS:-"8 16 32 64 128 256 512"}
EXTRA=${EXTRA:-""}             # explicit n:rows:cols runs appended to the grid
SCHEMES=${SCHEMES:-"2 3 1"}    # erase-program, vanilla, differential

# planeDecoder MUST sit immediately after sslSwitchMatrix -- it is the last
# entry of rw_main's P[] array, and the emit loop runs over P[] before
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

# ---------------------------------------------------------------- header ----
# UNCHANGED from sweep.sh.  Order must match rw_main.cpp's emit exactly; see
# the FIX 1 note there (nHalfSel is emitted LAST in the write group).
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
  printf ',eWr_cells_J,eWr_inhibit_J,eWr_lines_J'
  printf ',driveFrac,switchesPerWrite,writesToFailure,eWriteTrue_J,tWriteTrue_s'
  printf ',nHalfSel'
  # disturb margin
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
  printf '\n'
} > "${OUT}.csv"

# ------------------------------------------------------------- run list -----
# n:rows:cols triples -- the square grid, then the explicit extras
RUNS=""
for n in $NS; do
  for s in $SIZES; do RUNS="$RUNS ${n}:${s}:${s}"; done
done
RUNS="$RUNS $EXTRA"

n_cfg=0; for r in $RUNS; do n_cfg=$((n_cfg+1)); done
n_sch=0; for s in $SCHEMES; do n_sch=$((n_sch+1)); done
n_pat=0; for p in $PATTERNS; do n_pat=$((n_pat+1)); done
echo "square sweep: sizes [$SIZES] x planes [$NS]${EXTRA:+ + extra [$EXTRA]}"
echo "  $n_cfg configs x $n_sch schemes x $n_pat patterns = $((n_cfg*n_sch*n_pat)) runs -> ${OUT}.csv"

# ----------------------------------------------------------------- runs -----
n_runs=0; n_fail=0
for r in $RUNS; do
  N=${r%%:*}; tmp=${r#*:}; ROWS=${tmp%%:*}; COLS=${tmp##*:}
  for s in $SCHEMES; do
    for p in $PATTERNS; do
      P1=${p%%:*}; Q1=${p##*:}
      TAG="n${N}_${ROWS}x${COLS}_p${P1}_q${Q1}_s${s}"
      ./rw_main "$CELL" "$RDO" "$N" "$ROWS" "$COLS" "$P1" "$Q1" "$s" 0 "$VMIN" \
           > "${LOGS}/${TAG}.txt" 2>> "${OUT}.csv"
      rc=$?
      if [ $rc -gt 2 ]; then echo "  FAILED rc=$rc  $TAG"; n_fail=$((n_fail+1)); fi
      n_runs=$((n_runs+1))
    done
  done
  printf '  done n=%-4s %sx%s\n' "$N" "$ROWS" "$COLS"
done

NCOL=$(head -1 "${OUT}.csv" | awk -F, '{print NF}')
NROW=$(($(wc -l < "${OUT}.csv") - 1))
echo
echo "$n_runs runs ($n_fail hard failures), $NROW rows x $NCOL columns -> ${OUT}.csv   (full reports in ${LOGS}/)"
echo

# Column indices are looked up BY NAME so the tables below survive any future
# header change.
col() { head -1 "${OUT}.csv" | tr ',' '\n' | grep -n "^$1\$" | cut -d: -f1; }
I_RFJ=$(col eReadMacro_fJb);  I_WFJ=$(col eWriteMacro_fJb)
I_TR=$(col tReadMacro_s);     I_TW=$(col tWriteMacro_s)
I_AM=$(col areaMacro_um2);    I_AC=$(col areaCore_um2)
I_DV=$(col dVsense_V);        I_SV=$(col senseVerdict)
I_SW=$(col switchesPerWrite); I_WTF=$(col writesToFailure)
I_HS=$(col nHalfSel);         I_VH=$(col vHalfSel_V);  I_DR=$(col disturbRatio)

# ------------------------------------------------- grids: size down, n across
# One table per metric, vanilla scheme, random data.  These are the numbers
# the paper figures are drawn from (plot_rw_sweep.py reads the CSV directly).
grid() {   # $1 = column index, $2 = title, $3 = scale, $4 = printf fmt
  awk -F, -v c="$1" -v title="$2" -v k="$3" -v fmt="$4" -v ns="$NS" -v sizes="$SIZES" '
    NR>1 && $6=="0.5" && $8=="vanilla" && $4==$5 { v[$4","$3] = $c * k }
    END {
      nn = split(ns, N, " "); ns_ = split(sizes, S, " ")
      printf "\n  %s\n", title
      printf "  %8s |", "S \\ n"
      for (j=1;j<=nn;j++) printf " %9s", N[j]; printf "\n  %8s-+", "--------"
      for (j=1;j<=nn;j++) printf "----------"; printf "\n"
      for (i=1;i<=ns_;i++) {
        printf "  %8s |", S[i]
        for (j=1;j<=nn;j++) { key = S[i]","N[j]
          if (key in v) printf " " fmt, v[key]; else printf " %9s", "-" }
        printf "\n"
      }
    }' "${OUT}.csv"
}
echo "=== vanilla write, random data (p1 = q1 = 0.5); rows = cols = S ==="
grid "$I_RFJ" "read energy, macro  [fJ/bit]"        1     "%9.2f"
grid "$I_WFJ" "write energy, macro [fJ/bit]"        1     "%9.2f"
grid "$I_TR"  "read latency, macro  [ns]"           1e9   "%9.2f"
grid "$I_TW"  "write latency, macro [ns]"           1e9   "%9.2f"
grid "$I_AM"  "macro area  [um^2]  (negative = NeuroSim area model failed)" 1 "%9.0f"
grid "$I_DV"  "sense signal dV  [mV]  (VMIN = $VMIN V)" 1e3 "%9.2f"

# density and array efficiency are derived, so they get their own awk
awk -F, -v am="$I_AM" -v ac="$I_AC" -v ns="$NS" -v sizes="$SIZES" '
  NR>1 && $6=="0.5" && $8=="vanilla" && $4==$5 && $am>0 {
      d[$4","$3] = $3*$4*$5/$am; e[$4","$3] = 100*$ac/$am }
  END {
    nn = split(ns, N, " "); ns_ = split(sizes, S, " ")
    for (pass=1; pass<=2; pass++) {
      printf "\n  %s\n", (pass==1 ? "bit density  [Mbit/mm^2]  (all n planes)" : "array efficiency  [%]  = core / macro")
      printf "  %8s |", "S \\ n"; for (j=1;j<=nn;j++) printf " %9s", N[j]; printf "\n  %8s-+", "--------"
      for (j=1;j<=nn;j++) printf "----------"; printf "\n"
      for (i=1;i<=ns_;i++) { printf "  %8s |", S[i]
        for (j=1;j<=nn;j++) { key=S[i]","N[j]
          if (key in d) printf " %9.1f", (pass==1 ? d[key] : e[key]); else printf " %9s", "-" }
        printf "\n" } } }' "${OUT}.csv"

# ------------------------------------------------- sense verdict per n ------
# dV = Qsw/C_node is geometry-independent on a flat pillar, so one line per n.
echo
echo "  sense verdict vs n (vanilla, largest size):"
awk -F, -v dv="$I_DV" -v sv="$I_SV" 'NR>1 && $6=="0.5" && $8=="vanilla" && $4==$5 {
    if (!($3 in best) || $4 > best[$3]) { best[$3]=$4; v[$3]=$dv*1e3; s[$3]=$sv } }
  END { for (n in v) printf "  %6d  %8.2f mV  %s\n", n, v[n], s[n] }' "${OUT}.csv" | sort -n

# ------------------------------------------------- write-scheme comparison ---
echo
echo "=== write schemes, random data ==="
printf '%5s %6s | %23s | %17s | %s\n' \
       "n" "S" "     write fJ/bit      " " switches/write  " " writes to failure"
printf '%5s %6s | %7s %7s %7s | %5s %5s %5s | %9s %9s %9s\n' \
       "" "" "erase" "vanil" "diff" "erase" "vanil" "diff" "erase" "vanil" "diff"
awk -F, -v wf="$I_WFJ" -v sw="$I_SW" -v wtf="$I_WTF" 'NR>1 && $6=="0.5" && $4==$5 {
    k = sprintf("%05d %05d", $3, $4)
    wfb[k "|" $8] = $wf; s[k "|" $8] = $sw; w[k "|" $8] = $wtf; seen[k] = 1
}
END{
    for (k in seen) { split(k, f, " ")
        printf "%5d %6d | %7.2f %7.2f %7.2f | %5.3f %5.3f %5.3f | %9.2e %9.2e %9.2e\n",
               f[1], f[2],
               wfb[k"|erase-program"], wfb[k"|vanilla"], wfb[k"|differential"],
               s  [k"|erase-program"], s  [k"|vanilla"], s  [k"|differential"],
               w  [k"|erase-program"], w  [k"|vanilla"], w  [k"|differential"]
    }
}' "${OUT}.csv" | sort -k1,1n -k2,2n

# ------------------------------------------------- disturb margin -----------
# V_half/V_coercive > 1 means the V/2 scheme is WRITING the half-selected cells.
echo
echo "=== disturb margin (vanilla) ==="
awk -F, -v hs="$I_HS" -v vh="$I_VH" -v dr="$I_DR" '
NR==1 { printf "%5s %6s | %11s %9s %9s  %s\n", "n","S","half-sel","V_half","V_h/V_c","verdict"; next }
$6=="0.5" && $8=="vanilla" && $4==$5 {
    printf "%5d %6d | %11d %9.3f %9.3f  %s\n",
           $3,$4, $hs, $vh, $dr, ($dr > 1.0 ? "OVER COERCIVE -- not inhibiting" : "ok")
}' "${OUT}.csv" | sort -k1,1n -k2,2n -s

cat <<'NOTE'

  Next:  python3 build_rw_xlsx.py  <csv>            # workbook with Grid_* pivot sheets
         python3 plot_rw_sweep.py  <csv> --outdir figs   # paper figures (PDF + PNG)
NOTE
