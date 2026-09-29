/*******************************************************************************
 * mem_main.cpp  --  MEMORY OPERATOR MODE driver for DNN+NeuroSim V1.4 (3D FeRAM fork)
 *
 * WHAT THIS IS
 *   A second entry point that stands beside main.cpp.  It instantiates ONE
 *   SubArray, forces SINGLE-ROW ACTIVATION, and reports read/write/area PPA for
 *   the CORE ARRAY ONLY - cells and array lines, with no decoders, switch
 *   matrices, level shifters, sense amps, ADC, mux, shift-add or adder trees.
 *
 *   It deliberately does NOT go through Chip -> Tile -> ProcessingUnit, because
 *   that path is a dataflow mapper: it needs a netStructure file, it tiles a
 *   weight matrix, and it derives activityRowRead from an input vector.  A
 *   memory access has no weight matrix and no input vector.  Bolting a memory
 *   mode into ProcessingUnit would mean defeating the mapper at a dozen points;
 *   a separate main defeats it at zero.
 *
 * DEFINITION OF ONE ACCESS
 *   1 access = 1 word line, on 1 capacitor plane = numCol cells = numCol bits.
 *   With the default 512-column subarray that is a 512-bit page = 64 B = one
 *   cache line.  The n stacked planes are SEPARATE pages that share a pillar;
 *   they are not read together.  Every number this driver prints is per page.
 *
 * PREREQUISITES
 *   Apply, in order, the patch sets in MEMORY_MODE_GUIDE.md:
 *     Part 1  - eight correctness fixes (these change existing CIM numbers too)
 *     Part 2  - Param.h / Param.cpp fields
 *     Part 3  - SubArray.h core accumulators
 *     Part 4  - SubArray.cpp core/periphery split
 *
 * BUILD
 *   Put this file next to main.cpp.  In the makefile, build a SECOND binary:
 *       mem_main: $(filter-out main.o, $(OBJS)) mem_main.o
 *   Do not link main.o and mem_main.o together - both define main(), and both
 *   pull in Definition.h, which defines the global objects.
 *
 * USAGE
 *   ./mem_main <cellType> <readMode> <n> [numRow] [numCol]
 *     cellType : 1 SRAM | 2 RRAM | 5 2T-nC | 6 1T-nC | 7 1T1C   (ProcessingUnit.cpp:75-86)
 *     readMode : 1 ndro | 2 qndro | 3 dro                        (ProcessingUnit.cpp:95-101)
 *     n        : stacked capacitor planes; forced to 1 for planar cells
 *
 *   Human-readable table -> stdout.  One CSV line -> stderr, so a sweep can do
 *       ./mem_main 5 1 64 >/dev/null 2>>memcore.csv
 ******************************************************************************/
#include <cstdio>
#include <random>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <algorithm>

#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "SubArray.h"
#include "Definition.h"     /* defines the globals: inputParameter, tech, cell, param, gen */

using namespace std;

/* ---------------------------------------------------------------------------
 * Column resistance vector.
 * CalculateLatency/CalculatePower take one resistance per column.  For the
 * 2T-nC / 1T-nC branches this vector is only consumed by the sense-amp and
 * ADC models, which memory mode excludes, so a uniform ON-state column is the
 * neutral choice.  It DOES matter for the SRAM/RRAM baselines, where the
 * column current sets the sensing energy - there we want the ON state so the
 * baseline is not flattered.
 * ------------------------------------------------------------------------- */
static vector<double> ColumnResistance(int numCol, double rOn)
{
    /* For 2T-nC the sensed quantity is the RBL current through the read
     * transistor. The effective column resistance must be consistent with the
     * read model's own numbers: V_RBL / I_on = 0.5 V / 2 uA = 250 kohm.
     * Passing cell.resistanceOn (6 kohm) overstates sense energy by ~7x.   */
    const double V_RBL = 0.5, I_ON = 2e-6;      // SubArray.cpp 2T-nC read block
    return vector<double>(numCol, V_RBL / I_ON);
	
	//return vector<double>(numCol, rOn);
}

static void Row(const string &label, double v, const string &unit)
{
    cout << "  " << left << setw(32) << label
         << right << setw(16) << fixed << setprecision(4) << v
         << "  " << unit << '\n';
}

static void Head(const string &s)
{
    cout << '\n' << s << '\n' << string(70, '-') << '\n';
}

/* Emit one CSV field.  `valid` is false for cell types where the quantity
 * has no meaning (e.g. staircase area on a planar cell); the field is then
 * left empty so that every row still has the same number of columns. */
static void CsvField(double v, bool valid)
{
	if (valid) cerr << v;
	cerr << ',';
}

int main(int argc, char *argv[])
{
    /* ====================================================================
     * 1.  ARGUMENTS
     * ==================================================================== */
    const int cellType = (argc > 1) ? atoi(argv[1]) : 5;
    const int readMode = (argc > 2) ? atoi(argv[2]) : 1;
    int       nPlanes  = (argc > 3) ? atoi(argv[3]) : 64;

    const bool is3D = (cellType == 5 || cellType == 6);
    if (!is3D) nPlanes = 1;                 /* a planar cell has one plane */

    /* ====================================================================
     * 2.  PARAM
     *     Param's constructor has already run (Definition.h news it up), so
     *     every default from Param.cpp is in place.  We override only what a
     *     memory access changes.
     * ==================================================================== */
    param->memcelltype = cellType;
    param->mem_rdo     = readMode;
    switch (param->mem_rdo) {
        case 1: param->readVoltage = 0.1; break;   // ndro
        case 2: param->readVoltage = 1.2; break;   // qndro
        case 3: param->readVoltage = 2.0; break;   // dro
    }
    
    param->bitsPerCell = nPlanes;

    if (argc > 4) param->numRowSubArrayPhysical = atoi(argv[4]);
    if (argc > 5) param->numColSubArray         = atoi(argv[5]);

    /* numRowSubArray is the LOGICAL row count (pillars x planes).  Param.cpp:217
     * derives it from bitsPerCell, so it must be re-derived after we change n. */
    param->numRowSubArray = param->numRowSubArrayPhysical * param->bitsPerCell;

    /* ---- memory operator mode (Part 2 of the guide) ---- */
    param->memoryMode      = true;
    param->coreOnly        = true;
    param->numRowActivated = 1;             /* one word line per access */

    /* ---- switch off every CIM feature: a memory access has none of them ---- */
    param->operationmode      = 1;          /* conventionalSequential           */
    param->conventionalSequential = 1;
    param->conventionalParallel   = 0;
    param->BNNsequentialMode = param->BNNparallelMode  = 0;
    param->XNORsequentialMode = param->XNORparallelMode = 0;
    param->parallelRead       = 0;
    param->numRowParallel     = 1;
    param->numColMuxed        = 1;          /* no column multiplexing           */
    param->numBitInput        = 1;          /* one read pulse, not a DAC ramp   */
    param->synapseBit         = 1;
    param->cellBit            = 1;
    param->numColPerSynapse   = 1;
    param->numRowPerSynapse   = 1;
    param->SARADC             = false;
    param->levelOutput        = 2;          /* 1-bit sensing                    */
    param->synchronous        = false;      /* report seconds, not clock cycles */

    const int numRow = param->numRowSubArrayPhysical;
    const int numCol = param->numColSubArray;

    /* ====================================================================
     * 3.  TECHNOLOGY AND CELL
     *     Mirrors ProcessingUnitInitialize() (ProcessingUnit.cpp:75-200) so
     *     the cell this driver builds is bit-identical to the CIM one.
     * ==================================================================== */
    switch (param->memcelltype) {
        case 7: cell.memCellType = Type::_1T1C; break;
        case 6: cell.memCellType = Type::_1TnC; break;
        case 5: cell.memCellType = Type::_2TnC; break;
        case 4: cell.memCellType = Type::Cap;   break;
        case 3: cell.memCellType = Type::FeFET; break;
        case 2: cell.memCellType = Type::RRAM;  break;
        case 1: cell.memCellType = Type::SRAM;  break;
        default: cerr << "bad cellType\n"; return 1;
    }
    switch (param->accesstype) {
        case 4: cell.accessType = none_access;  break;
        case 3: cell.accessType = diode_access; break;
        case 2: cell.accessType = BJT_access;   break;
        case 1: cell.accessType = CMOS_access;  break;
        default: break;
    }
    switch (param->mem_rdo) {
        case 3: cell.mem_rdo = Type::dro;   break;
        case 2: cell.mem_rdo = Type::qndro; break;
        case 1: cell.mem_rdo = Type::ndro;  break;
        default: break;
    }
    switch (param->transistortype) {
        case 3: inputParameter.transistorType = TFET;         break;
        case 2: inputParameter.transistorType = FET_2D;       break;
        default: inputParameter.transistorType = conventional; break;
    }
    inputParameter.deviceRoadmap = (param->deviceroadmap == 1) ? HP : LSTP;
    inputParameter.temperature   = param->temp;
    inputParameter.processNode   = param->technode;
    tech.Initialize(inputParameter.processNode,
                    inputParameter.deviceRoadmap,
                    inputParameter.transistorType);

    cell.resistanceOn     = param->resistanceOn;
    cell.resistanceOff    = param->resistanceOff;
    cell.resistanceAvg    = (cell.resistanceOn + cell.resistanceOff) / 2;
    cell.readVoltage      = param->readVoltage;
    cell.readPulseWidth   = param->readPulseWidth;
    cell.accessVoltage    = param->accessVoltage;
    cell.resistanceAccess = param->resistanceAccess;
    cell.featureSize      = param->featuresize;
    cell.writeVoltage     = param->writeVoltage;
    cell.writePulseWidth  = param->writePulseWidth;

    const char *cn = "";
    if (cell.memCellType == Type::SRAM) {
        cn = "SRAM 6T";
        cell.heightInFeatureSize = param->heightInFeatureSizeSRAM;
        cell.widthInFeatureSize  = param->widthInFeatureSizeSRAM;
        cell.widthSRAMCellNMOS   = param->widthSRAMCellNMOS;
        cell.widthSRAMCellPMOS   = param->widthSRAMCellPMOS;
        cell.widthAccessCMOS     = param->widthAccessCMOS;
        cell.minSenseVoltage     = param->minSenseVoltage;
    } else if (cell.memCellType == Type::_2TnC) {
        cn = "2T-nC";
        cell.heightInFeatureSize = param->heightInFeatureSize2TnC;
        cell.widthInFeatureSize  = param->widthInFeatureSize2TnC;
    } else if (cell.memCellType == Type::_1TnC) {
        cn = "1T-nC";
        cell.heightInFeatureSize = param->heightInFeatureSize1TnC;
        cell.widthInFeatureSize  = param->widthInFeatureSize1TnC;
    } else if (cell.memCellType == Type::_1T1C) {
        cn = "1T1C planar FeRAM";
        cell.heightInFeatureSize = param->heightInFeatureSize1T1C;
        cell.widthInFeatureSize  = param->widthInFeatureSize1T1C;
    } else {
        cn = "1T1R (RRAM)";
        cell.heightInFeatureSize = (cell.accessType == CMOS_access)
                                 ? param->heightInFeatureSize1T1R
                                 : param->heightInFeatureSizeCrossbar;
        cell.widthInFeatureSize  = (cell.accessType == CMOS_access)
                                 ? param->widthInFeatureSize1T1R
                                 : param->widthInFeatureSizeCrossbar;
    }

    const char *rn = (!is3D && cell.memCellType != Type::_1T1C) ? "n/a"
                   : (readMode == 1 ? "ndro" : readMode == 2 ? "qndro" : "dro");

    /* ====================================================================
     * 4.  SUBARRAY
     * ==================================================================== */
    SubArray *sub = new SubArray(inputParameter, tech, cell);

    sub->conventionalSequential = true;
    sub->conventionalParallel   = false;
    sub->BNNsequentialMode  = sub->BNNparallelMode  = false;
    sub->XNORsequentialMode = sub->XNORparallelMode = false;
    sub->SARADC      = false;
    sub->levelOutput = 2;
    sub->validated   = false;   /* no 1.4x beta fudge on a plain memory access */
    sub->numColMuxed       = 1;
    sub->numReadPulse      = 1;
    sub->numCellPerSynapse = 1;
    sub->avgWeightBit      = 1;
    sub->maxNumWritePulse  = 1;
    sub->numWritePulse     = 1;
    sub->numRowParallel    = 1;
    sub->numAdd            = 1;
    sub->relaxArrayCellHeight = param->relaxArrayCellHeight;
    sub->relaxArrayCellWidth  = param->relaxArrayCellWidth;
    sub->currentMode       = param->currentMode;
    sub->spikingMode       = NONSPIKING;

    /* ------------------------------------------------------------------
     * SINGLE-ROW ACTIVATION.  This is the heart of memory mode.
     *
     * Every array expression in SubArray.cpp is written in the shape
     *     <per-cell term> * numRow * activityRow{Read,Write} * numCol
     * i.e. "per-cell cost x number of activated rows x row width".  Setting
     * activityRow* = 1/numRow makes numRow*activityRow* == 1 exactly, so the
     * SAME formula now describes one asserted word line.  No formula is
     * rewritten; the activation factor is what changes.
     *
     * activityRowWrite is READ at ~30 sites in SubArray.cpp and ASSIGNED at
     * none, in SubArray.cpp or ProcessingUnit.cpp.  It is an uninitialised
     * double.  Setting it here is not a convenience, it is a bug fix.
     * ------------------------------------------------------------------ */
    sub->activityRowRead  = 1.0 / numRow;
    sub->activityRowWrite = 1.0 / numRow;
    sub->activityColRead  = 1.0;   /* the whole page is sensed  */
    sub->activityColWrite = 1.0;   /* the whole page is written */

    sub->numReadCellPerOperationNeuro   = numCol;
    sub->numWriteCellPerOperationNeuro  = numCol;
    sub->numReadCellPerOperationMemory  = numCol;
    sub->numWriteCellPerOperationMemory = numCol;
    sub->numReadCellPerOperationFPGA    = numCol;
    sub->numWriteCellPerOperationFPGA   = numCol;

    sub->clkFreq = param->clkFreq;        // must be set BEFORE Initialize()
    sub->Initialize(numRow, numCol, param->unitLengthWireResistance);
    sub->CalculateArea();

    vector<double> colRes = ColumnResistance(numCol, cell.resistanceOn);

    /* ------------------------------------------------------------------
     * Two-pass clock, exactly as main.cpp:299-316 does it.  Pass 1 with
     * CalculateclkFreq=true measures the critical path; pass 2 uses it.
     * With param->synchronous == false the second pass reports seconds, which
     * is what we want - a memory macro's latency is a time, not a cycle count.
     *
     * writeBack=false: in memory mode the restore is NOT an amortised refresh
     * decided by a loop counter in ProcessingUnit.  It is a property of the
     * read-out mode, and Part 4 of the guide moves it inside SubArray where it
     * belongs.  Passing true here would double-count it.
     * ------------------------------------------------------------------ */
    sub->CalculateLatency(1e20, colRes, true,  false);
    if (sub->readLatency > 0) param->clkFreq = 1.0 / sub->readLatency;
    sub->clkFreq = param->clkFreq;
    sub->CalculateLatency(1e20, colRes, false, false);
    sub->CalculatePower(colRes, false);

    /* ====================================================================
     * 5.  EXTRACTION
     * ==================================================================== */
    double eRead, eWrite, tRead, tWrite, areaCore;

    if (is3D) {
        eRead    = sub->readEnergyCore;
        eWrite   = sub->writeEnergyCore;
        tRead    = sub->readLatencyCore;
        tWrite   = sub->writeLatencyCore;
        areaCore = sub->areaCore;
    } else {
        /* SRAM / RRAM / 1T1C: the stock *Array variables are already array-only
         * (they are never contaminated with write-back periphery), so read them
         * directly.  There is no core/periphery latency split for these - see
         * the note printed below. */
        eRead    = sub->readDynamicEnergyArray;
        eWrite   = sub->writeDynamicEnergyArray;
        tRead    = sub->readLatency;
        tWrite   = sub->writeLatency;
        areaCore = sub->areaArray;
    }

    const double bits     = numCol;
    const double capacity = (double)numRow * numCol * nPlanes;   /* bits */

    cout << "\n======================================================================\n"
         << "  " << cn << "   read-out: " << rn << "   n = " << nPlanes
         << "   " << param->technode << " nm\n"
         << "  MEMORY OPERATOR MODE - one word line per access, core array only\n"
         << "======================================================================\n";

    Head("ORGANISATION");
    cout << "  " << left << setw(32) << "subarray"
         << right << setw(16) << (to_string(numRow) + " x " + to_string(numCol)) << '\n';
    Row("planes (n)",                nPlanes, "");
    Row("rows asserted per access",  numRow * sub->activityRowRead, "(must be 1)");
    Row("page width",                bits, "bits");
    Row("pages in subarray",         (double)numRow * nPlanes, "");
    Row("capacity",                  capacity / 8.0 / 1048576.0, "MB");

    Head("AREA");
    Row("array",        sub->areaArray * 1e12, "um^2");
    if (is3D) Row("staircase contacts", sub->areaStaircase * 1e12, "um^2");
    Row("CORE TOTAL",   areaCore * 1e12, "um^2");
    Row("per Mbit",     areaCore * 1e12 / (capacity / 1048576.0), "um^2/Mbit");
    Row("bit density",  capacity / (areaCore * 1e12) * 1e6 / 1e9, "Gb/mm^2");

    Head("LATENCY, CORE ONLY");
    if (is3D) {
        Row("charge transfer",  param->chargeDelay * 1e9, "ns");
        Row("line setup RC",    sub->lineSetupCore * 1e9, "ns");
        Row("restore",          sub->restoreLatencyCore * 1e9, "ns");
        Row("READ  core",       tRead  * 1e9, "ns");
        Row("WRITE core",       tWrite * 1e9, "ns");
        Row("read  bandwidth",  bits / tRead  / 1e9, "Gb/s");
        Row("write bandwidth",  bits / tWrite / 1e9, "Gb/s");
        Head("LATENCY, WITH PERIPHERY (reference)");
        Row("READ  total",  sub->readLatency  * 1e9, "ns");
        Row("WRITE total",  sub->writeLatency * 1e9, "ns");
    } else {
        Row("READ  (incl. periphery)",  tRead  * 1e9, "ns");
        Row("WRITE (incl. periphery)",  tWrite * 1e9, "ns");
        cout << "  NOTE: this cell type has no core/periphery latency split.\n"
             << "        Do not put these on the same axis as a 3D core latency.\n";
    }

    Head("ENERGY, CORE ONLY");
    if (is3D) {
        Row("read  sense",      sub->readSenseEnergyCore   * 1e12, "pJ");
        Row("read  restore",    sub->readRestoreEnergyCore * 1e12, "pJ");
        Row("write cells",      sub->writeCellEnergyCore   * 1e12, "pJ");
        Row("write inhibition", sub->inhibitionEnergyCore  * 1e12, "pJ");
    }
    Row("READ  per access",  eRead  * 1e12, "pJ");
    Row("WRITE per access",  eWrite * 1e12, "pJ");
    Row("read  per bit",     eRead  / bits * 1e15, "fJ/bit");
    Row("write per bit",     eWrite / bits * 1e15, "fJ/bit");

    if (is3D) {
        const double restorePeriph = sub->readDynamicEnergyArray - sub->readEnergyCore;
        Head("READ ENERGY — WHAT EACH NUMBER CONTAINS");
        Row("core array (cells + lines)",        sub->readEnergyCore * 1e12, "pJ");
        Row("+ restore periphery",               restorePeriph * 1e12, "pJ");
        Row("= array + restore",                 sub->readDynamicEnergyArray * 1e12, "pJ");
        cout << "  Restore periphery = the WBL/WWL/WPL/SSL decoder and switch-matrix\n"
             << "  energy of the write-back that a DESTRUCTIVE read triggers. It is zero\n"
             << "  for ndro, which performs no restore, so for ndro the two rows above are\n"
             << "  equal by construction. Neither figure includes the row drivers, level\n"
             << "  shifters, sense amplifier, mux or ADC of a normal read.\n";
    }

    const double eSense = sub->senseEnergyCore;
    const double tSense = sub->senseLatencyCore;

    Head("LATENCY, SENSE TIER  (core array + column sense amp, no row drivers)");
    Row("core array",      tRead * 1e9, "ns");
    Row("sense amplifier", tSense * 1e9, "ns");
    Row("READ core+sense", (tRead + tSense) * 1e9, "ns");

    Head("ENERGY, SENSE TIER");
    Row("core array",      eRead * 1e12, "pJ");
    Row("sense amplifier", eSense * 1e12, "pJ");
    Row("READ core+sense", (eRead + eSense) * 1e12, "pJ");
    Row("per bit",         (eRead + eSense) / bits * 1e15, "fJ/bit");

    /* ---- guard rails ------------------------------------------------- */
    bool bad = false;
    if (eRead <= 0 || eWrite <= 0) {
        cout << "\n  *** a core energy is <= 0: the corresponding formula is\n"
                "  *** still commented out in SubArray.cpp.  DO NOT PLOT.\n";
        bad = true;
    }
    if (is3D && fabs(numRow * sub->activityRowRead - 1.0) > 1e-9) {
        cout << "\n  *** activation is not single-row.  DO NOT PLOT.\n";
        bad = true;
    }
    //if (is3D && tRead >= sub->readLatency && sub->readLatency > 0) {
    if (is3D && (tRead - sub->restoreLatencyCore) >= sub->readLatency && sub->readLatency > 0) {  
        cout << "\n  *** core read latency is not smaller than the full read\n"
                "  *** latency: the Part 4 latency patch did not take effect.\n";
        bad = true;
    }
    cout << (bad ? "\n  STATUS: FAILED\n" : "\n  STATUS: ok\n") << endl;

    if (is3D && sub->senseEnergyCore <= 0) {
    	cout << "\n  *** senseEnergyCore is 0: the sense-tier patch is missing from\n"
    	        "  *** the " << cn << " branch of CalculatePower().  The core+sense\n"
    	        "  *** columns are meaningless.  DO NOT PLOT.\n";
    	bad = true;
    }
    if (is3D && sub->writeCellEnergyCore <= 0) {
    	cout << "\n  *** writeCellEnergyCore is 0: the charge-based write patch is\n"
    	        "  *** missing from the " << cn << " branch.  DO NOT PLOT.\n";
    	bad = true;
    }

    	/* 1T-nC / 1T1C read destructively dumps 2Pr*A onto the bit line.  Q/C_BL is
	 * an UPPER BOUND on the signal (the ferroelectric cap is in the denominator
	 * too).  If even this bound is under minSenseVoltage, the latency and energy
	 * columns describe a read that never resolves.  2T-nC is a gain cell read in
	 * the current domain, so the quantity does not apply and is left blank. */
	const bool isChargeShare = (cell.memCellType == Type::_1TnC ||
	                            cell.memCellType == Type::_1T1C);
	double vSignal = 0.0;
	if (isChargeShare && sub->capBL > 0)
		vSignal = (param->twoPr * param->cellAreaFE) / sub->capBL;

    	
    // /* CSV -> stderr */
    // cerr << cn << ',' << rn << ',' << nPlanes << ','
    //      << numRow << ',' << numCol << ','
    //      << areaCore * 1e12 << ','
    //      << tRead * 1e9 << ',' << tWrite * 1e9 << ','
    //      << eRead * 1e12 << ',' << eWrite * 1e12 << ','
    //      << eRead / bits * 1e15 << ',' << eWrite / bits * 1e15 << '\n';

    // delete sub;
    // return bad ? 1 : 0;
    //

    	/* ====================================================================
	 * 6.  CSV LINE  ->  stderr
	 * One line per run, no header.  The header is written once by the sweep
	 * script so that many runs can be appended to a single file.
	 * ==================================================================== */
	{
		/* Derived here, not in SubArray.cpp, so that nothing depends on the
		 * order in which CalculatePower() assigns its members. */
		double eRestorePeriph = is3D ? (sub->readDynamicEnergyArray - eRead) : 0.0;
		double eArrayPlusRest = is3D ?  sub->readDynamicEnergyArray          : eRead;
		double tSenseTier     = is3D ? (tRead + sub->senseLatencyCore)       : 0.0;
		double eSenseTier     = is3D ? (eRead + sub->senseEnergyCore)        : 0.0;

		cerr << setprecision(8);

		cerr << cn << ',' << rn << ',' << nPlanes << ','
		     << numRow << ',' << numCol << ',';
		CsvField(capacity / 1048576.0,                     true);   /* Mbit */

		/* ---- area ---- */
		CsvField(sub->areaArray     * 1e12,                true);
		CsvField(sub->areaStaircase * 1e12,                is3D);
		CsvField(areaCore           * 1e12,                true);
		CsvField(areaCore * 1e12 / (capacity / 1048576.0), true);
		CsvField(capacity / (areaCore * 1e12) * 1e6 / 1e9, true);

		/* ---- latency, core only ---- */
		CsvField(param->chargeDelay      * 1e9,            is3D);
		CsvField(sub->lineSetupCore      * 1e9,            is3D);
		CsvField(sub->restoreLatencyCore * 1e9,            is3D);
		CsvField(tRead                   * 1e9,            true);
		CsvField(tWrite                  * 1e9,            true);
		CsvField(bits / tRead  / 1e9,                      true);
		CsvField(bits / tWrite / 1e9,                      true);

		/* ---- latency, full macro (reference) ---- */
		CsvField(sub->readLatency  * 1e9,                  true);
		CsvField(sub->writeLatency * 1e9,                  true);

		/* ---- energy, core only ---- */
		CsvField(sub->readSenseEnergyCore   * 1e12,        is3D);
		CsvField(sub->readRestoreEnergyCore * 1e12,        is3D);
		CsvField(sub->writeCellEnergyCore   * 1e12,        is3D);
		CsvField(sub->inhibitionEnergyCore  * 1e12,        is3D);
		CsvField(eRead  * 1e12,                            true);
		CsvField(eWrite * 1e12,                            true);
		CsvField(eRead  / bits * 1e15,                     true);
		CsvField(eWrite / bits * 1e15,                     true);

		/* ---- what the read number contains ---- */
		// CsvField(eRestorePeriph * 1e12,                    is3D);
		// CsvField(eArrayPlusRest * 1e12,                    true);

		CsvField(sub->readDynamicEnergyArray * 1e12, true);   /* full-macro array figure, as the model computes it */

		/* ---- sense tier ---- */
		CsvField(sub->senseLatencyCore * 1e9,              is3D);
		CsvField(tSenseTier            * 1e9,              is3D);
		CsvField(sub->senseEnergyCore  * 1e12,             is3D);
		CsvField(eSenseTier            * 1e12,             is3D);
		CsvField(eSenseTier / bits * 1e15,                 is3D);

		CsvField(vSignal * 1e3,              isChargeShare);
		CsvField(param->minSenseVoltage*1e3, isChargeShare);

		cerr << (bad ? "FAILED" : "ok") << '\n';
	}
}
