#!/usr/bin/env bash
# sweep_new.sh -- FULL FACTORIAL GRID, VANILLA WRITE SCHEME ONLY.
#
#   planes n                 : 64, 256, 512
#   rows                     : 64, 128, 512
#   cols                     : 64, 128, 512
#   write scheme             : 3 = vanilla single plate cycle  (ONLY)
#
#   3 x 3 x 3 = 27 runs.  Three distinct values on every axis, so every trend
#   can be fitted AND its curvature tested -- two points can only ever draw a
#   straight line.
#
# WHY THESE VALUES
#   cols  64/128/512   128 sits on the fitted area optimum (~120 columns) that
#                      the aspect-ratio model predicted by extrapolating between
#                      the two corners.  This is the run that confirms or kills
#                      it.  It also makes tPlateRC ~ numCol^2 a three-point
#                      demonstration: 7.15 / 28.59 / 457.5 ps, exactly x4, x16.
#   rows  64/128/512   with cols 64/128/512 this yields THREE iso-capacity
#                      points at 32768 cells -- 64x512, 128x256 (not in grid),
#                      512x64 -- and 128x128 / 512x512 as the square anchors.
#                      NOTE: the row direction is nearly flat over this span
#                      (2.6% -> 16.6% far-row current loss).  The row WALL is
#                      past ~1000 rows; run that separately with
#                      ./sweep.sh 5 2 rows, which scans 64 -> 8192.
#   n     64/256/512   256 -> 512 crosses the sense wall.  dV = Qsw/capSense is
#                      29.37 mV at n=256 and 14.70 mV at n=512 against a 25 mV
#                      vMinSense, so n=512 is the first FAILING height and the
#                      pair brackets the limit at n ~ 300.  n=512 is also where
#                      the plate driver passes 50% of the read.
#
# Single-row read and single-row write throughout: rw_main sets
# activityRowRead = activityRowWrite = 1/numRow, so numRow x activity == 1.
#
# Usage:
#   ./sweep_new.sh             # 2T-nC qndro, random data, vanilla only
#   ./sweep_new.sh 6 3         # 1T-nC dro
#   ./sweep_new.sh 5 2 corners # + the all-0 / all-1 data corners (x5 runs)
#
# Output:  rw_<cell>_<rdo>_grid.csv   ->  python3 build_grid_xlsx.py

set -u
CELL=${1:-5}
RDO=${2:-2}
MODE=${3:-plain}
VMIN=0.025                     # offset-cancelled sense amp; 0 = use Param's value
OUT=rw_${CELL}_${RDO}_grid
LOGS=${OUT}_logs
mkdir -p "$LOGS"

[ -x ./rw_main ] || { echo "build it first:  make rw_main"; exit 1; }

# ------------------------------------------------------------ the grid ------
NS="64 256 512"
ROWS_SET="64 256 512"
COLS_SET="64 256 512"
EXTRA=""                       # explicit n:rows:cols runs appended to the grid.
                               # 8:8:8 is deliberately NOT here: below ~16 rows
                               # NeuroSim's area model fails (four row-direction
                               # decoders each return ~-9.8e8 um^2) and the run
                               # reports a negative macro area.

# FIX 4: planeDecoder MUST sit immediately after sslSwitchMatrix -- it is the
# last entry of rw_main's P[] array, and the emit loop runs over P[] before
# appending the sense amp / latch / cellArray groups by hand.  The same order
# is hardcoded in build_grid_xlsx.py's BLOCKS list.
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
SCHEMES="3"                    # VANILLA ONLY -- this is the whole point of this
                               # script.  The three-scheme comparison lives in
                               # sweep.sh; repeating it here would triple the
                               # grid for a result that does not depend on shape
                               # (the V/2 inhibit is shape-blind).

# ---------------------------------------------------------------- header ----
# IDENTICAL to sweep.sh.  The order must match rw_main.cpp's emit EXACTLY --
# csv.DictReader zips names to POSITIONS, so a single inserted name rotates
# every field after it.  Do not reorder anything here; append only.
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
  # LAST in this group, not fourth.
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
RUNS=""
for n in $NS; do
  for r in $ROWS_SET; do
    for c in $COLS_SET; do
      RUNS="$RUNS ${n}:${r}:${c}"
    done
  done
done
RUNS="$RUNS $EXTRA"

# ----------------------------------------------------------------- runs -----
n_runs=0
for r in $RUNS; do
  [ -z "$r" ] && continue
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
NDATA=$(( $(wc -l < "${OUT}.csv") - 1 ))
echo
echo "$n_runs runs, $NDATA data rows, $NCOL columns -> ${OUT}.csv   (reports in ${LOGS}/)"
[ "$NCOL" -ne 194 ] && echo "  WARNING: expected 194 columns.  Header and rw_main's emit are out of step."
echo

# ----------------------------------------------------- column index lookup ---
# Every table below resolves its fields BY NAME, so appending new diagnostics to
# the header never silently rotates a table again.
idx () { head -1 "${OUT}.csv" | tr ',' '\n' | grep -n "^$1\$" | cut -d: -f1; }
I_N=$(idx n);              I_R=$(idx rows);            I_C=$(idx cols)
I_ERD=$(idx eReadMacro_fJb);   I_EWR=$(idx eWriteMacro_fJb)
I_TRD=$(idx tReadMacro_s);     I_TWR=$(idx tWriteMacro_s)
I_TRC=$(idx tReadCore_s);      I_TWC=$(idx tWriteCore_s)
I_AC=$(idx areaCore_um2);      I_AM=$(idx areaMacro_um2);  I_LK=$(idx leakMacro_W)
I_DV=$(idx dVsense_V);         I_SV=$(idx senseVerdict)
I_INH=$(idx eWr_inhibit_J);    I_CN=$(idx capNode_F)
I_TPL=$(idx t_plateRC_s);      I_TSC=$(idx t_settleCol_s); I_TSR=$(idx t_settleRow_s)
I_TSD=$(idx t_senseDev_s);     I_TSA=$(idx t_senseAmp_s)
I_CC=$(idx capColSense_F);     I_DL=$(idx degenLoss_frac); I_VR=$(idx rblDropRatio)
I_DR=$(idx disturbRatio);      I_HS=$(idx nHalfSel);       I_VH=$(idx vHalfSel_V)

# --------------------------------------------------------- headline table ---
echo "HEADLINE  (vanilla, p1=q1=0.5)"
awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v erd=$I_ERD -v ewr=$I_EWR -v trd=$I_TRD \
        -v twr=$I_TWR -v am=$I_AM -v dv=$I_DV -v sv=$I_SV '
BEGIN{printf "%5s %6s %6s  %10s %12s  %9s %10s  %11s %8s %7s\n",
      "n","rows","cols","read fJ/b","write fJ/b","tRead ns","tWrite ns","macro um2","dV mV","sense"}
NR>1 && $6=="0.5" && $8=="vanilla" {
    printf "%5d %6d %6d  %10.3f %12.3f  %9.3f %10.3f  %11.1f %8.2f %7s\n",
           $n,$r,$c, $erd, $ewr, $trd*1e9, $twr*1e9, $am, $dv*1e3, $sv
}' "${OUT}.csv" | { IFS= read -r hdr; echo "$hdr"; sort -k1,1n -k2,2n -k3,3n; }

# ------------------------------------------------- GRID: rows x cols per n ---
# The pivot the whole sweep exists for.  One block per n, rows down, cols across.
# POSIX awk only -- no asorti/gensym.  The axis values come in from the shell,
# already in order, so nothing needs sorting inside awk.
grid () {   # $1 = field index, $2 = scale, $3 = printf fmt, $4 = title
  echo
  echo "$4"
  awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v f=$1 -v sc=$2 -v fmt="$3" \
          -v NSET="$NS" -v RSET="$ROWS_SET" -v CSET="$COLS_SET" '
  BEGIN{ nn=split(NSET,ns," "); nr=split(RSET,rs," "); nc=split(CSET,cs," ") }
  NR>1 && $6=="0.5" && $8=="vanilla" { v[$n" "$r" "$c]=$f*sc; seen[$n" "$r" "$c]=1 }
  END{
      for (i=1;i<=nn;i++) {
          printf "  n = %d\n      rows\\cols", ns[i]
          for (k=1;k<=nc;k++) printf "%12d", cs[k]
          printf "%14s\n", "x across cols"
          for (j=1;j<=nr;j++) {
              printf "  %10d", rs[j]
              first=0; last=0
              for (k=1;k<=nc;k++) {
                  key = ns[i]" "rs[j]" "cs[k]
                  if (!(key in seen)) { printf "%12s", "-"; continue }
                  x = v[key]
                  printf "%12"fmt, x
                  if (k==1) first=x
                  last=x
              }
              printf "%13.2fx\n", (first!=0 ? last/first : 0)
          }
          printf "  %10s", "x down rows"
          for (k=1;k<=nc;k++) {
              a=v[ns[i]" "rs[1]" "cs[k]]; b=v[ns[i]" "rs[nr]" "cs[k]]
              printf "%11.2fx", (a!=0 ? b/a : 0)
          }
          printf "\n\n"
      }
  }' "${OUT}.csv"
}
grid $I_TRD 1e9   ".4f" "GRID  read latency, macro (ns)         [row trend = across cols; col trend = down rows]"
grid $I_ERD 1     ".3f" "GRID  read energy (fJ/bit)"
grid $I_EWR 1     ".2f" "GRID  write energy (fJ/bit)"
grid $I_AM  1     ".1f" "GRID  macro area (um^2)"

# ------------------------------------------------------ iso-capacity check ---
# 64x512, 128x128 and 512x64 are not iso-capacity, but 64x512 and 512x64 are
# (32768 cells) and 128x128 (16384) is the square anchor.  The aspect-ratio
# model predicted the AREA minimum near 120 columns; if that is real, 512x64
# must beat 64x512 and the gap must be ~1.74x.
echo "ISO-CAPACITY  32768 cells: 64x512 (wide) vs 512x64 (tall)"
awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v am=$I_AM -v erd=$I_ERD -v ewr=$I_EWR \
        -v trd=$I_TRD -v NSET="$NS" '
BEGIN{ nn=split(NSET,ns," ") }
NR>1 && $6=="0.5" && $8=="vanilla" && (($r==64&&$c==512)||($r==512&&$c==64)) {
    k=$n; if ($r==64) { aw[k]=$am; ew[k]=$erd; ww[k]=$ewr; tw[k]=$trd }
          else        { at[k]=$am; et[k]=$erd; wt[k]=$ewr; tt[k]=$trd }
}
END{
  printf "%5s | %11s %11s %9s | %9s %9s %7s | %10s %10s %7s\n",
         "n","wide um2","tall um2","tall/wide","wide fJ/b","tall fJ/b","ratio","wide wfJ/b","tall wfJ/b","ratio"
  for (i=1;i<=nn;i++){ k=ns[i]
    if (!(k in aw) || !(k in at)) continue
    printf "%5d | %11.1f %11.1f %9.3f | %9.3f %9.3f %7.3f | %10.2f %10.2f %7.3f\n",
      k, aw[k], at[k], at[k]/aw[k], ew[k], et[k], et[k]/ew[k], ww[k], wt[k], wt[k]/ww[k]
  }
}' "${OUT}.csv"

# ------------------------------------------------------- n trend (planes) ----
# capSense = n*C_FE + C_gate, so dV ~ 1/n and the V/2 inhibit ~ (n-1).  Both are
# shape-blind, so any one geometry shows the law; 128x128 is used here.
echo
echo "PLANE TREND  at 128x128  (dV = Qsw/capNode; inhibit per bit is shape-blind)"
awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v cn=$I_CN -v dv=$I_DV -v sv=$I_SV \
        -v inh=$I_INH -v trd=$I_TRD -v ewr=$I_EWR '
BEGIN{printf "%6s | %10s %9s %10s | %11s %11s %8s\n",
      "n","capNode fF","dV mV","dV x n","inhib fJ/b","write fJ/b","sense"}
NR>1 && $6=="0.5" && $8=="vanilla" && $r==128 && $c==128 {
    printf "%6d | %10.2f %9.2f %10.2f | %11.2f %11.2f %8s\n",
      $n, $cn*1e15, $dv*1e3, $dv*1e3*$n, $inh*1e15/128, $ewr, $sv
}' "${OUT}.csv" | { IFS= read -r h; echo "$h"; sort -k1,1n; }
cat <<'NOTE'
  dV x n should be flat: dV = Qsw/(n*C_FE + C_gate) is 1/n to within the fixed
  C_gate.  inhibit/bit should follow (n-1)*1.2016 fJ at EVERY geometry.
NOTE

# ------------------------------------------------- disturb margin (FIX 7) ---
echo
awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v hs=$I_HS -v vh=$I_VH -v dr=$I_DR '
BEGIN{printf "%5s %6s %6s | %11s %9s %9s  %s\n",
      "n","rows","cols","half-sel","V_half","V_h/V_c","verdict"}
NR>1 && $6=="0.5" && $8=="vanilla" {
    printf "%5d %6d %6d | %11d %9.3f %9.3f  %s\n",
      $n,$r,$c, $hs, $vh, $dr, ($dr > 1.0 ? "OVER COERCIVE -- not inhibiting" : "ok")
}' "${OUT}.csv" | { IFS= read -r h; echo "$h"; sort -k1,1n -k2,2n -k3,3n; }

# ----------------------------------------------------------- column path -----
# The row count enters the read through the COLUMN line (the read transistor's
# source line): its capacitance sets the column settle + development time and
# its resistance degenerates the far row.  rblDropRatio > 1 means the far
# columns lose V_DS and a page-wide current read is not viable at that width.
echo
echo "COLUMN PATH  (vanilla)"
awk -F, -v n=$I_N -v r=$I_R -v c=$I_C -v cc=$I_CC -v tsc=$I_TSC -v tsd=$I_TSD \
        -v tsa=$I_TSA -v dl=$I_DL -v vr=$I_VR -v trd=$I_TRD -v tpl=$I_TPL '
BEGIN{printf "%5s %6s %6s | %8s %10s %10s %10s | %9s %9s | %10s  %s\n",
      "n","rows","cols","C_col fF","t_plate ns","t_col ns","t_sense ns","far-row","RBL drop","tRead ns","verdict"}
NR>1 && $6=="0.5" && $8=="vanilla" {
    printf "%5d %6d %6d | %8.2f %10.4f %10.4f %10.4f | %8.2f%% %9.3f | %10.3f  %s\n",
      $n,$r,$c, $cc*1e15, $tpl*1e9, $tsc*1e9, $tsa*1e9, 100*$dl, $vr, $trd*1e9,
      ($vr > 1.0 ? "RBL IR: page current collapses V_DS" : "ok")
}' "${OUT}.csv" | { IFS= read -r h; echo "$h"; sort -k1,1n -k2,2n -k3,3n; }

echo
echo "next:  python3 build_grid_xlsx.py ${OUT}.csv"
