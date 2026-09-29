#!/bin/bash
# 2T-nC, all three read-out modes, n = 8..256, full field set.
OUT=sweep_2TnC_new.csv
mkdir -p logs

cat > "$OUT" <<'HDR'
cell,mode,n,numRow,numCol,capacity_Mbit,areaArray_um2,areaStaircase_um2,areaCore_um2,area_um2_per_Mbit,density_Gb_per_mm2,tCharge_ns,tLineSetup_ns,tRestore_ns,tReadCore_ns,tWriteCore_ns,readBW_Gbps,writeBW_Gbps,tReadFull_ns,tWriteFull_ns,eReadSense_pJ,eReadRestore_pJ,eWriteCells_pJ,eWriteInhibit_pJ,eReadCore_pJ,eWriteCore_pJ,eReadCore_fJ_per_bit,eWriteCore_fJ_per_bit,eRestorePeriph_pJ,eArrayPlusRestore_pJ,tSenseAmp_ns,tReadCoreSense_ns,eSenseAmp_pJ,eReadCoreSense_pJ,eReadCoreSense_fJ_per_bit,status
HDR

for m in 1 2 3; do                          # 1=ndro 2=qndro 3=dro
  for n in 8 16 32 64 128 256; do
    ./mem_main 5 "$m" "$n" > "logs/2TnC_m${m}_n${n}.txt" 2>> "$OUT"
  done
done

# column-count check: header and every data row must agree
awk -F, 'NR==1{h=NF} NF!=h{print "COLUMN MISMATCH on line "NR": "NF" vs "h; e=1} END{if(!e) print "columns ok ("h")"}' "$OUT"
grep -c FAILED "$OUT" | xargs -I{} echo "rows with STATUS=FAILED: {}"
echo "wrote $OUT"
