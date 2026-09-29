			/* ---- block-level accounting (see BlockStats.h) ---------------------
			 * One hook, at the end of the 2T-nC power branch.  ProcessingUnit
			 * calls CalculateLatency() immediately before CalculatePower() on
			 * this same object, so every .readLatency / .writeLatency read below
			 * is current, and this only runs on the !CalculateclkFreq pass.
			 * Energy sums reproduce the chip total exactly; latency sums are
			 * total block-time, NOT the MAX-reduced critical path.            */
			{
			    BlockStats &B = BLK();
			    B.enabled = true;
			    B.calls  += 1;
			    B.bitsPerCell = BlockStats::mx(B.bitsPerCell, (double)param->bitsPerCell);
#define X(n) B.rdLat_##n += n.readLatency;  B.rdEn_##n += n.readDynamicEnergy;  \
             B.wrLat_##n += n.writeLatency; B.wrEn_##n += n.writeDynamicEnergy; \
             B.leak_##n   = BlockStats::mx(B.leak_##n, n.leakage);
			    BLK_COMPONENTS(X)
#undef X
#define X(n) B.tot_##n += (double)(n);
			    BLK_TOTALS(X)
#undef X
#define X(n) B.core_##n += (double)(n);
			    BLK_CORE(X)
#undef X
#define X(n) B.cfg_##n = BlockStats::mx(B.cfg_##n, (double)(n));
			    BLK_CONFIG(X)
#undef X
			    B.maxReadLatency  = BlockStats::mx(B.maxReadLatency,  readLatency);
			    B.maxWriteLatency = BlockStats::mx(B.maxWriteLatency, writeLatency);
			}
			/* -------------------- end block-level accounting -------------------- */
