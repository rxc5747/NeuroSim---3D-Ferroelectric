/*******************************************************************************
 * BlockStats.h -- chip-level, per-block latency / energy accounting
 *
 * WHY THIS AND NOT cout IN SubArray
 *   A cout inside SubArray::CalculatePower fires once per (tile, PE, subarray,
 *   input-vector) -- thousands of lines per layer, and no single line is the
 *   answer.  The chip total is built by a chain of += :
 *       Chip.cpp:798        *readDynamicEnergy += tileReadDynamicEnergy;
 *       Tile.cpp            *readDynamicEnergy += peReadDynamicEnergy;
 *       ProcessingUnit:378  *readDynamicEnergy += subArray->readDynamicEnergy;
 *   and every level loops over its instances, calling the level below once per
 *   instance.  So summing a per-block value at the SAME point where
 *   ProcessingUnit adds subArray->readDynamicEnergy reproduces the chip total
 *   exactly, with the same multiplicity.  That is what this does.
 *
 * WHERE IT HOOKS
 *   One call, at the end of the 2T-nC branch of SubArray::CalculatePower.
 *   ProcessingUnit.cpp calls CalculateLatency() immediately before
 *   CalculatePower() on the same object, so every component's .readLatency and
 *   .writeLatency are still current at that point -- one hook covers both.
 *   It also only runs on the !CalculateclkFreq pass, i.e. pass 2, which is the
 *   pass whose numbers are real.  Pass 1 only measures the critical path.
 *
 * HOW TO READ THE OUTPUT -- the one caveat that matters
 *   ENERGY is summed the same way the chip sums it, so the per-block energies
 *   ADD UP to the chip read/write energy.  Trust them.
 *   LATENCY is NOT additive at chip level: Tile and Chip combine child
 *   latencies with MAX(), not +, because tiles run in parallel.  The latency
 *   numbers here are therefore TOTAL BLOCK-TIME summed over every invocation --
 *   the right thing for "where does the array spend its time", and NOT equal to
 *   the chip's reported readLatency.  Both are printed so you can see the gap.
 *   LEAKAGE is a power, so it is tracked as a MAX (the per-subarray value),
 *   not a sum.
 *
 * Printing happens automatically at program exit (static destructor), so there
 * is nothing to add to main.cpp.
 ******************************************************************************/
#ifndef BLOCKSTATS_H_
#define BLOCKSTATS_H_

#include <cstdio>
#include <cstdlib>
#include <string>

/* Every block we want a line for.  Add or remove a name here and the struct,
 * the accumulator and the printout all follow automatically. */
#define BLK_COMPONENTS(X)                                                      \
    X(rslDecoder)  X(rslSwitchMatrix)  X(rblDecoder)  X(rblSwitchMatrix)       \
    X(wblDecoder)  X(wblSwitchMatrix)  X(wwlDecoder)  X(wwlSwitchMatrix)       \
    X(wplDecoder)  X(wplSwitchMatrix)  X(sslDecoder)  X(sslSwitchMatrix)       \
    X(mux)         X(muxDecoder)                                               \
    X(sarADC)      X(multilevelSenseAmp) X(multilevelSAEncoder)                \
    X(currentSenseAmp)                                                         \
    X(shiftAddWeight) X(shiftAddInput)  X(adder)      X(dff)

/* SubArray-level roll-ups (summed) */
#define BLK_TOTALS(X)                                                          \
    X(readLatency) X(readLatencyADC) X(readLatencyAccum) X(readLatencyOther)   \
    X(writeLatency)                                                            \
    X(readDynamicEnergy) X(readDynamicEnergyArray)                             \
    X(readDynamicEnergyADC) X(readDynamicEnergyAccum) X(readDynamicEnergyOther)\
    X(writeDynamicEnergy) X(writeDynamicEnergyArray)

/* Core-tier members added by the 3D FeRAM fork (summed).  Delete a line here
 * if your SubArray.h does not have that member. */
#define BLK_CORE(X)                                                            \
    X(readEnergyCore) X(readSenseEnergyCore) X(readRestoreEnergyCore)          \
    X(senseEnergyCore)                                                         \
    X(writeEnergyCore) X(writeCellEnergyCore) X(inhibitionEnergyCore)          \
    X(writeLatencyCore) X(restoreLatencyCore) X(lineSetupCore)

/* Configuration constants -- recorded as MAX, they do not vary per call */
#define BLK_CONFIG(X)                                                          \
    X(capWBL) X(capWWL) X(capWPL) X(capSSL) X(capRSL) X(capRBL)                \
    X(resRow) X(resCol) X(numRow) X(numCol)

struct BlockStats {
#define X(n) double rdLat_##n, rdEn_##n, wrLat_##n, wrEn_##n, leak_##n;
    BLK_COMPONENTS(X)
#undef X
#define X(n) double tot_##n;
    BLK_TOTALS(X)
#undef X
#define X(n) double core_##n;
    BLK_CORE(X)
#undef X
#define X(n) double cfg_##n;
    BLK_CONFIG(X)
#undef X
    double calls;
    double maxReadLatency, maxWriteLatency;
    double bitsPerCell;
    std::string tag;
    bool enabled;

    BlockStats() { reset(); enabled = false; tag = ""; }

    void reset() {
#define X(n) rdLat_##n = rdEn_##n = wrLat_##n = wrEn_##n = leak_##n = 0;
        BLK_COMPONENTS(X)
#undef X
#define X(n) tot_##n = 0;
        BLK_TOTALS(X)
#undef X
#define X(n) core_##n = 0;
        BLK_CORE(X)
#undef X
#define X(n) cfg_##n = 0;
        BLK_CONFIG(X)
#undef X
        calls = 0; maxReadLatency = 0; maxWriteLatency = 0; bitsPerCell = 0;
    }

    static double mx(double a, double b) { return a > b ? a : b; }

    ~BlockStats() { print(); }

    void print() const {
        if (!enabled || calls == 0) return;
        const double eR = tot_readDynamicEnergy;
        const double eW = tot_writeDynamicEnergy;
        double lR = 0;
#define X(n) lR += rdLat_##n;
        BLK_COMPONENTS(X)
#undef X

        printf("\n");
        printf("==========================================================================\n");
        const char *envTag = getenv("BLOCKSTATS_TAG");
        printf(" BLOCK-LEVEL TOTALS %s\n", envTag ? envTag : tag.c_str());
        printf(" summed over %.0f subarray invocations (whole network, pass 2 only)\n", calls);
        printf("==========================================================================\n");
        printf(" config: numRow=%.0f numCol=%.0f bitsPerCell=%.0f\n",
               cfg_numRow, cfg_numCol, bitsPerCell);
        printf("         capWBL=%.4g capWWL=%.4g capWPL=%.4g capSSL=%.4g capRSL=%.4g capRBL=%.4g F\n",
               cfg_capWBL, cfg_capWWL, cfg_capWPL, cfg_capSSL, cfg_capRSL, cfg_capRBL);
        printf("         resRow=%.4g resCol=%.4g ohm\n", cfg_resRow, cfg_resCol);

        printf("\n--- READ ------------------------------------------------------------------\n");
        printf(" %-22s %14s %7s %14s %7s\n", "block", "energy [pJ]", "%", "time [ns]", "%");
#define X(n) if (rdEn_##n != 0 || rdLat_##n != 0)                              \
        printf(" %-22s %14.6g %6.2f%% %14.6g %6.2f%%\n", #n,                   \
               rdEn_##n * 1e12,  (eR > 0 ? 100.0 * rdEn_##n / eR  : 0.0),      \
               rdLat_##n * 1e9,  (lR > 0 ? 100.0 * rdLat_##n / lR : 0.0));
        BLK_COMPONENTS(X)
#undef X
        printf(" %-22s %14.6g %6.2f%% %14.6g %6.2f%%\n", "[cell array]",
               tot_readDynamicEnergyArray * 1e12,
               (eR > 0 ? 100.0 * tot_readDynamicEnergyArray / eR : 0.0), 0.0, 0.0);
        printf(" %-22s %14s %7s %14s\n", "--- sum of blocks", "", "", "");
        printf(" %-22s %14.6g %6s  %14.6g\n", "subarray total", eR * 1e12, "100%", lR * 1e9);
        printf(" %-22s %14s %7s %14.6g   <- MAX, not the sum\n",
               "chip critical path", "", "", maxReadLatency * 1e9);

        printf("\n--- WRITE -----------------------------------------------------------------\n");
        printf(" %-22s %14s %7s %14s\n", "block", "energy [pJ]", "%", "time [ns]");
#define X(n) if (wrEn_##n != 0 || wrLat_##n != 0)                              \
        printf(" %-22s %14.6g %6.2f%% %14.6g\n", #n,                           \
               wrEn_##n * 1e12, (eW > 0 ? 100.0 * wrEn_##n / eW : 0.0),        \
               wrLat_##n * 1e9);
        BLK_COMPONENTS(X)
#undef X
        printf(" %-22s %14.6g %6.2f%%\n", "[cell array]",
               tot_writeDynamicEnergyArray * 1e12,
               (eW > 0 ? 100.0 * tot_writeDynamicEnergyArray / eW : 0.0));
        printf(" %-22s %14.6g %6s  %14.6g\n", "subarray total", eW * 1e12, "100%",
               tot_writeLatency * 1e9);

        printf("\n--- SUBARRAY ROLL-UPS (the numbers NeuroSim itself reports) ---------------\n");
#define X(n) printf(" %-34s %16.6g\n", #n, tot_##n);
        BLK_TOTALS(X)
#undef X

        printf("\n--- CORE TIER (3D FeRAM fork members) -------------------------------------\n");
#define X(n) printf(" %-34s %16.6g\n", #n, core_##n);
        BLK_CORE(X)
#undef X

        printf("\n--- LEAKAGE, per subarray (W) ---------------------------------------------\n");
#define X(n) if (leak_##n != 0) printf(" %-22s %16.6g\n", #n, leak_##n);
        BLK_COMPONENTS(X)
#undef X
        printf("==========================================================================\n\n");
        fflush(stdout);
    }
};

/* Single instance.  Its destructor runs at program exit, so the report prints
 * itself -- no main.cpp edit needed. */
inline BlockStats& BLK() { static BlockStats s; return s; }

#endif  /* BLOCKSTATS_H_ */
