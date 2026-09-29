#!/bin/bash
# 2T-nC and 1T-nC, all three read-out modes, n = 8..256.
OUT=sweep_mem_new_1.csv
mkdir -p logs

cat > "$OUT" <<'HDR'
cell,mode,n,numRow,numCol,capacity_Mbit,areaArray_um2,areaStaircase_um2,areaCore_um2,area_um2_per_Mbit,density_Gb_per_mm2,tCharge_ns,tLineSetup_ns,tRestore_ns,tReadCore_ns,tWriteCore_ns,readBW_Gbps,writeBW_Gbps,tReadFull_ns,tWriteFull_ns,eReadSense_pJ,eReadRestore_pJ,eWriteCells_pJ,eWriteInhibit_pJ,eReadCore_pJ,eWriteCore_pJ,eReadCore_fJ_per_bit,eWriteCore_fJ_per_bit,eReadArrayFull_pJ,tSenseAmp_ns,tReadCoreSense_ns,eSenseAmp_pJ,eReadCoreSense_pJ,eReadCoreSense_fJ_per_bit,vSignal_mV,minSense_mV,status
HDR

for c in 5 6; do                            # 5 = 2T-nC, 6 = 1T-nC
  for m in 2 3; do                        # 1 = ndro, 2 = qndro, 3 = dro
    for n in 64 256; do
      ./mem_main "$c" "$m" "$n" > "logs/c${c}_m${m}_n${n}.txt" 2>> "$OUT"
    done
  done
done

awk -F, 'NR==1{h=NF} NF!=h{print "COLUMN MISMATCH line "NR": "NF" vs "h; e=1}
         END{if(!e) print "columns ok ("h")"}' "$OUT"
echo "FAILED rows: $(grep -c ',FAILED$' "$OUT")"
grep ',FAILED$' "$OUT" | cut -d, -f1-3
echo "wrote $OUT"
