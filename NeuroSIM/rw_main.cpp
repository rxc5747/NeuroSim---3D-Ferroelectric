/*******************************************************************************
 * rw_main.cpp - READ / WRITE OPERATOR MODE 
 *
 *       line activation -> charge transfer -> bit-line settle -> sense -> latch
 *
 *   INCLUDED periphery:
 *       row / plane address decoders      rblDecoder, rslDecoder, wblDecoder,
 *                                         wwlDecoder, wplDecoder, sslDecoder
 *       line drivers (input drivers)      rblSwitchMatrix, rslSwitchMatrix,
 *                                         wblSwitchMatrix, wwlSwitchMatrix,
 *                                         wplSwitchMatrix, sslSwitchMatrix
 *       readout periphery                 currentSenseAmp (1 per column)
 *                                         output latch (DFF, 1 per I/O bit)
 *
 *   EXCLUDED periphery (CIM-only):
 *       sarADC, multilevelSenseAmp, multilevelSAEncoder,
 *       mux + muxDecoder (column sharing for the ADC),
 *       shiftAddWeight, shiftAddInput, adder, accumulation dff,
 *       the beta = 1.4 CIM validation multiplier.
 *
 * USAGE
 *   ./rw_main <cell> <rdo> <n> <rows> <cols> [p1] [q1] [scheme] [ioWidth]
 *     cell   : 5 = 2T-nC | 6 = 1T-nC | 7 = 1T1C | 1 = SRAM | 2 = RRAM
 *     rdo    : 1 = ndro | 2 = qndro | 3 = dro
 *     n      : stacked capacitor planes (1 for planar cells)
 *     rows   : PHYSICAL pillar rows in the subarray
 *     cols   : columns  (= page width in bits)
 *     p1     : fraction of the READ page holding the switching state ('1').
 *              default 0.5.  .
 *     q1     : fraction of the WRITTEN word that is '1'.  default 0.5.
 *     ioWidth: bits driven off the macro per access.  default = cols (page out).
 *
 ******************************************************************************/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>          
#include <sstream>
#include <string>
#include <vector>

#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "DFF.h"
#include "SubArray.h"
#include "RowDecoder.h"   
#include "Definition.h"     

using namespace std;

static const double EPS0_VACUUM = 8.854e-12;  // F/m, vacuum permittivity 
static const double V_RBL_BIAS  = 0.5;        // RBL/RSL DC bias during read 

static const double R_PLATE_DRV = 5.0e3;      

static const double T_WRITE_RRAM = 50e-9;     

static double DegenerateCurrent(double vov, double rPath, double iTh, double ss, double iOn)
{
    double I = MIN(iTh * pow(10.0, vov / ss), iOn);
    if (rPath <= 0.0) return I;
    for (int k = 0; k < 200; k++) {
        const double In    = MIN(iTh * pow(10.0, (vov - I * rPath) / ss), iOn);
        const double Inext = 0.5 * (I + In);
        if (fabs(Inext - I) < 1e-9 * MAX(I, 1e-15)) { I = Inext; break; }
        I = Inext;
    }
    return I;
}

static void H(const string &s)
{ cout << '\n' << s << '\n' << string(78, '-') << '\n'; }

static void R(const string &label, double v, const string &unit, int prec = 4)
{
    cout << "  " << left << setw(38) << label
         << right << setw(16) << scientific << setprecision(prec) << v
         << "  " << unit << '\n';
}

static void Rf(const string &label, double v, const string &unit, int prec = 4)
{
    cout << "  " << left << setw(38) << label
         << right << setw(16) << fixed << setprecision(prec) << v
         << "  " << unit << '\n';
}

static void Note(const string &s) { cout << "    " << s << '\n'; }

int main(int argc, char *argv[])
{
    //  1.  ARGUMENTS
    const int cellType = (argc >  1) ? atoi(argv[1]) : 5;
    const int readMode = (argc >  2) ? atoi(argv[2]) : 2;
    int       nPlanes  = (argc >  3) ? atoi(argv[3]) : 64;
    const int argRows  = (argc >  4) ? atoi(argv[4]) : 512;
    const int argCols  = (argc >  5) ? atoi(argv[5]) : 512;
    const double p1    = (argc >  6) ? atof(argv[6]) : 0.5;
    const double q1    = (argc >  7) ? atof(argv[7]) : 0.5;
    const int scheme   = (argc >  8) ? atoi(argv[8]) : 2;
    const int argIO    = (argc >  9) ? atoi(argv[9]) : 0;
    const double argVmin = (argc > 10) ? atof(argv[10]) : 0.0;
    const double argVread = (argc > 11) ? atof(argv[11]) : 0.0;

    if (p1 < 0 || p1 > 1 || q1 < 0 || q1 > 1) {
        cerr << "p1 and q1 must be in [0,1]\n"; return 1;
    }
    if (scheme < 1 || scheme > 3) {
        cerr << "scheme must be 1 (differential RMW), 2 (erase-then-program), "
                "or 3 (vanilla single plate cycle)\n";
        return 1;
    }

    const bool is3D = (cellType == 5 || cellType == 6);
    const bool isFeRAM = (cellType == 5 || cellType == 6 || cellType == 7);
    if (!is3D) nPlanes = 1;

    // 2.  PARAM - override only what a memory access changes
    param->memcelltype = cellType;
    param->mem_rdo     = readMode;
    param->bitsPerCell = nPlanes;
    param->numRowSubArrayPhysical = argRows;
    param->numColSubArray         = argCols;
    param->numRowSubArray = param->numRowSubArrayPhysical * (int)param->bitsPerCell;

    param->operationmode          = 1;      
    param->conventionalSequential = 1;
    param->conventionalParallel   = 0;
    param->BNNsequentialMode  = param->BNNparallelMode  = 0;
    param->XNORsequentialMode = param->XNORparallelMode = 0;
    param->parallelRead     = 0;
    param->numRowParallel   = 1;
    param->numColMuxed      = 1;            
    param->numBitInput      = 1;            
    param->synapseBit       = 1;
    param->cellBit          = 1;
    param->numColPerSynapse = 1;
    param->numRowPerSynapse = 1;
    param->SARADC           = false;
    param->levelOutput      = 2;            
    param->synchronous      = false;       

    const int numRow  = param->numRowSubArrayPhysical;
    const int numCol  = param->numColSubArray;
    const int ioWidth = (argIO > 0 && argIO <= numCol) ? argIO : numCol;

    // 3.  TECHNOLOGY AND CELL - mirrors ProcessingUnitInitialize()
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
        case 3: inputParameter.transistorType = TFET;          break;
        case 2: inputParameter.transistorType = FET_2D;        break;
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

    static const double R_ON_RRAM    = 100e3;   
    static const double R_ONOFF_RRAM = 17.0;    
    static const double V_ACCESS_RRAM = 1.1;    
    if (cellType == 2) {
        cell.resistanceOn   = R_ON_RRAM;
        cell.resistanceOff  = R_ON_RRAM * R_ONOFF_RRAM;
        cell.resistanceAvg  = 0.5 * (cell.resistanceOn + cell.resistanceOff);
        cell.resMemCellOn   = cell.resistanceOn  * (1.0 + IR_DROP_TOLERANCE);
        cell.resMemCellOff  = cell.resistanceOff + cell.resistanceOn * IR_DROP_TOLERANCE;
        cell.resMemCellAvg  = 0.5 * (cell.resMemCellOn + cell.resMemCellOff);
        cell.resCellAccess  = cell.resistanceOn * IR_DROP_TOLERANCE;
    }
    cell.readVoltage      = param->readVoltage;

    const bool   twoTpre    = (cell.memCellType == Type::_2TnC);
    const double capFEpre   = EPS0_VACUUM * param->epsFE
                            * param->cellAreaFE / param->tFE;
    const double capGatePre = capFEpre / param->capGateRatio;
    const double qswPre     = param->twoPr * param->cellAreaFE;
    const double dVmodPre   = qswPre / ((double)nPlanes * capFEpre + capGatePre);

    const double readWindow = param->ssReadTr
                            * log10(param->ionSatReadTr / param->ioffReadTr);
    const double readOv     = (param->readOverdriveAuto && twoTpre)
                            ? MAX(0.0, readWindow - 0.5 * dVmodPre)
                            : param->readOverdrive;

    {
        const double vCoercive = param->ecFE * 1e8 * param->tFE;   
        const double vBiasNode = param->vthReadTr + readOv;
        if      (cell.mem_rdo == Type::qndro) cell.readVoltage = 0.75*vCoercive + vBiasNode;
        else if (cell.mem_rdo == Type::ndro ) cell.readVoltage = 0.25*vCoercive + vBiasNode;
        else                                  cell.readVoltage = param->writeVoltage;
    }

    if (argVread > 0) cell.readVoltage = argVread;
    cell.readPulseWidth   = param->readPulseWidth;
    cell.accessVoltage    = (cellType == 2) ? V_ACCESS_RRAM   
                                           : param->accessVoltage;
    cell.resistanceAccess = param->resistanceAccess;
    cell.featureSize      = param->featuresize;
    cell.writeVoltage     = param->writeVoltage;
    cell.writePulseWidth  = param->writePulseWidth;
    if (isFeRAM) cell.widthAccessCMOS = param->widthAccessCMOS;

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
        cn = "2T-nC vertical 3D FeRAM";
        cell.heightInFeatureSize = param->heightInFeatureSize2TnC;
        cell.widthInFeatureSize  = param->widthInFeatureSize2TnC;
    } else if (cell.memCellType == Type::_1TnC) {
        cn = "1T-nC vertical 3D FeRAM";
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
    const char *rn = readMode == 1 ? "ndro" : readMode == 2 ? "qndro" : "dro";
    const char *sn = (scheme == 2) ? "erase-then-program (blind, 2 phases)"
                   : (scheme == 3) ? "vanilla single plate cycle (blind, 2 edges)"
                                   : "differential read-modify-write (owes a read)";

    // 4.  SUBARRAY - single-row activation
    SubArray *sub = new SubArray(inputParameter, tech, cell);

    sub->conventionalSequential = true;
    sub->conventionalParallel   = false;
    sub->BNNsequentialMode  = sub->BNNparallelMode  = false;
    sub->XNORsequentialMode = sub->XNORparallelMode = false;
    sub->SARADC      = false;
    sub->levelOutput = 2;
    sub->validated   = false;          
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

    // THE single-row activation.  numRow * activityRow* == 1 exactly. 
    sub->activityRowRead  = 1.0 / numRow;
    sub->activityRowWrite = 1.0 / numRow;
    sub->activityColRead  = 1.0;
    sub->activityColWrite = 1.0;

    sub->numReadCellPerOperationNeuro   = numCol;
    sub->numWriteCellPerOperationNeuro  = numCol;
    sub->numReadCellPerOperationMemory  = numCol;
    sub->numWriteCellPerOperationMemory = numCol;
    sub->numReadCellPerOperationFPGA    = numCol;
    sub->numWriteCellPerOperationFPGA   = numCol;

    if (param->clkFreq <= 0) param->clkFreq = 1e9;
    sub->clkFreq = param->clkFreq;

    sub->Initialize(numRow, numCol, param->unitLengthWireResistance);
    sub->CalculateArea();

    vector<double> colRes(numCol, cell.resistanceOn);

    sub->CalculateLatency(1e20, colRes, true,  false);
    if (sub->readLatency > 0) param->clkFreq = 1.0 / sub->readLatency;
    sub->clkFreq = param->clkFreq;
    sub->CalculateLatency(1e20, colRes, false, false);
    sub->CalculatePower(colRes, false);

    // 5.  LINE PARASITICS (read back from the SubArray)
    const double capWBL = sub->capWBL;    
    const double capRBL = sub->capRBL;    
    const double capRSL = sub->capRSL;    
    const double capWWL = sub->capWWL;    
    const double capWPL = sub->capWPL;    
    const double capSSL = sub->capSSL;    
    const double resRow = sub->resRow;
    const double resCol = sub->resCol;
    const double capWBLpar  = (sub->capWBLpar  > 0.0) ? sub->capWBLpar  : capWBL;
    const double capWBLwire = (sub->capWBLwire > 0.0) ? sub->capWBLwire : 0.0;

    const double rPlateMeasured = sub->wblSwitchMatrix.resTg;
    const bool   rPlateIsReal   = (rPlateMeasured > 1.0 && rPlateMeasured < 1.0e7);
    const double rPlateDrv      = rPlateIsReal ? rPlateMeasured : R_PLATE_DRV;

    const bool twoT = (cell.memCellType == Type::_2TnC);

    const double CAP_FE_LIN = EPS0_VACUUM * param->epsFE
                            * param->cellAreaFE / param->tFE;
    const double capGate = CAP_FE_LIN / param->capGateRatio;

    const double capSense = twoT
        ? ((double)nPlanes * CAP_FE_LIN + capGate)   
        : (capWBL);                                  
    const char *senseNodeName = twoT ? "pillar node (n*C_FE + parasitic)"
                                     : "bit line (capWBL)";
    const double capDriven = twoT ? capWBL : capWBL;

     // 6.  CELL PHYSICS - 1 versus 0 model
    const double Qsw = param->twoPr * param->cellAreaFE;       
    const double Vr  = cell.readVoltage;
    const double Vw  = cell.writeVoltage;

    const double vIntBias = param->vthReadTr + readOv;   
    const double vUnsel   = (twoT && nPlanes > 1)
                          ? (vIntBias * capSense - CAP_FE_LIN * Vr)
                            / ((nPlanes - 1) * CAP_FE_LIN)
                          : 0.0;
    const double dVunclamped = Qsw / capSense;   
    const double dVheadroom  = MAX(Vr - vIntBias, 0.0);
    const double dVmod       = MIN(dVunclamped, dVheadroom);
    const bool   dVclipped   = (dVunclamped > dVheadroom);

    const double iRead1 = MIN(param->ioffReadTr
                              * pow(10.0, (vIntBias + 0.5*dVmod - param->vthReadTr)
                                          / param->ssReadTr),
                              param->ionSatReadTr);
    const double iRead0 = MIN(param->ioffReadTr
                              * pow(10.0, (vIntBias - 0.5*dVmod - param->vthReadTr)
                                          / param->ssReadTr),
                              param->ionSatReadTr);

    double alphaRd = 0.0;                    
    double restoreDuty = 0.0;                
    if (cell.mem_rdo == Type::dro)   { alphaRd = 1.0; restoreDuty = 1.0; }
    else if (cell.mem_rdo == Type::qndro) {
        alphaRd     = param->qndroSwitchFraction;
        restoreDuty = 1.0 / param->qndroRefreshInterval;
    }

    const double eCellRead1 = alphaRd * Qsw * Vr + CAP_FE_LIN * Vr * Vr;
    const double eCellRead0 =                      CAP_FE_LIN * Vr * Vr;
    double eCellReadPage = numCol * (p1 * eCellRead1 + (1.0 - p1) * eCellRead0);

    double iDCPage = 0.0;
    if (twoT) iDCPage = numCol * (p1 * iRead1 + (1.0 - p1) * iRead0);

    const double nUnselStrips   = (twoT && nPlanes > 1 && param->readDrivesUnselPlanes)
                                ? ((double)nPlanes - 1.0) : 0.0;
    const double eReadUnselStrips = nUnselStrips * capWBLwire * vUnsel * vUnsel;
    double eReadLines = (twoT ? capWBLpar : capDriven) * Vr * Vr + eReadUnselStrips;
    if (twoT) eReadLines += capRBL * V_RBL_BIAS * V_RBL_BIAS
                          + capRSL * V_RBL_BIAS * V_RBL_BIAS;

    const double ePrecharge = 0.0;

    const double T_WRITE = param->writePulseWidth;   

    const double eWriteField = (Vw / param->tFE) / 1e8;              
    const double tauNLS      = param->nlsTauInf
                             * exp(pow(param->nlsEa / eWriteField, param->nlsAlpha));

    const double eSwitchCell = Qsw * Vw;
    const double eLinearCell = CAP_FE_LIN * Vw * Vw;
    const double eFlipCell   = eSwitchCell + eLinearCell;  
    
        const double vInh = Vw / param->inhibitDivider;               

    const double vCoercive2  = param->ecFE * 1e8 * param->tFE;      
    const double eHalfMVcm   = (vInh / param->tFE) / 1e8;           
    const double disturbRatio= (vCoercive2 > 0) ? vInh / vCoercive2 : 0.0;

    double nHalfSel = (double)numCol * ((double)nPlanes - 1.0);
    if (param->inhibitPillarModel == 1)
        nHalfSel += (double)numCol * ((double)numRow - 1.0) * (double)nPlanes;

    const double kInh     = param->inhibitSupplyDraw ? 1.0 : 0.5;
    const double eInhibit = kInh * CAP_FE_LIN * vInh * vInh * nHalfSel;
    
    
    const double eRestorePage = restoreDuty
                              * (p1 * numCol * eFlipCell + eInhibit);

    const double hamming = p1*(1.0-q1) + q1*(1.0-p1);

    double flipFrac;     
    double driveFrac;    
    double nWritePulses;
    switch (scheme) {
        case 2:  flipFrac = p1 + q1;  driveFrac = 2.0;      nWritePulses = 2; break;
        case 3:  flipFrac = hamming;  driveFrac = 1.0;      nWritePulses = 2; break;
        default: flipFrac = hamming;  driveFrac = hamming;  nWritePulses = 2; break;
    }
    const bool owesARead = (scheme == 1);

    const double eWriteCells   = numCol * (flipFrac  * eSwitchCell
                                         + driveFrac * eLinearCell);
    const double eWriteInhibit = nWritePulses * eInhibit;
    const double nUnselStripsW   = (twoT && nPlanes > 1) ? ((double)nPlanes - 1.0) : 0.0;
    const double eWriteUnselStrips = nWritePulses * nUnselStripsW * capWBLwire * vInh * vInh;
    const double eWriteLines   = nWritePulses *
                                 ((twoT ? capWBLpar : capWBL) + capWWL + capWPL + capSSL) * Vw * Vw
                               + eWriteUnselStrips;

     // 7.  PERIPHERY, read straight back out of the SubArray objects.
    const double WBL_ACT_FIX = 1.0 / (double)nPlanes;

    const double WBL_HW_FIX = (double)numCol / ((double)numRow * (double)nPlanes);

    const bool STRIP_SCAN_DFF = true;
#define SM_R(sm) ((sm).readDynamicEnergy  - (STRIP_SCAN_DFF ? (sm).dff.readDynamicEnergy : 0.0))
#define SM_W(sm) ((sm).writeDynamicEnergy - (STRIP_SCAN_DFF ? (sm).dff.readDynamicEnergy : 0.0))

    const bool havePlaneDecoder = (isFeRAM && nPlanes > 1);
    RowDecoder planeDecoder(inputParameter, tech, cell);
    if (havePlaneDecoder) {
        planeDecoder.Initialize(REGULAR_ROW,
                                (int)MAX(1.0, ceil(log2((double)nPlanes))),
                                false, false);
        planeDecoder.CalculateArea(NULL, NULL, NONE);
        planeDecoder.CalculateLatency(1e20, capWPL, NULL, resRow, numCol, 1, 1);
        planeDecoder.CalculatePower(1, 1);
    }

    struct Blk { const char *name; double eR, eW, tR, tW, lk, ar; };
    Blk P[] = {
      {"rblDecoder",      sub->rblDecoder.readDynamicEnergy,      sub->rblDecoder.writeDynamicEnergy,
                          sub->rblDecoder.readLatency,            sub->rblDecoder.writeLatency,
                          sub->rblDecoder.leakage,                sub->rblDecoder.area},
      {"rblSwitchMatrix", SM_R(sub->rblSwitchMatrix), SM_W(sub->rblSwitchMatrix),
                          sub->rblSwitchMatrix.readLatency,       sub->rblSwitchMatrix.writeLatency,
                          sub->rblSwitchMatrix.leakage,           sub->rblSwitchMatrix.area},
      {"rslDecoder",      sub->rslDecoder.readDynamicEnergy,      sub->rslDecoder.writeDynamicEnergy,
                          sub->rslDecoder.readLatency,            sub->rslDecoder.writeLatency,
                          sub->rslDecoder.leakage,                sub->rslDecoder.area},
      {"rslSwitchMatrix", SM_R(sub->rslSwitchMatrix), SM_W(sub->rslSwitchMatrix),
                          sub->rslSwitchMatrix.readLatency,       sub->rslSwitchMatrix.writeLatency,
                          sub->rslSwitchMatrix.leakage,           sub->rslSwitchMatrix.area},
      {"wblDecoder",      sub->wblDecoder.readDynamicEnergy,      sub->wblDecoder.writeDynamicEnergy,
                          sub->wblDecoder.readLatency,            sub->wblDecoder.writeLatency,
                          sub->wblDecoder.leakage,                sub->wblDecoder.area},
      {"wblSwitchMatrix", SM_R(sub->wblSwitchMatrix) * WBL_ACT_FIX,
                          SM_W(sub->wblSwitchMatrix) * WBL_ACT_FIX,
                          sub->wblSwitchMatrix.readLatency,       sub->wblSwitchMatrix.writeLatency,
                          sub->wblSwitchMatrix.leakage * WBL_HW_FIX,
                          sub->wblSwitchMatrix.area    * WBL_HW_FIX},
      {"wwlDecoder",      sub->wwlDecoder.readDynamicEnergy,      sub->wwlDecoder.writeDynamicEnergy,
                          sub->wwlDecoder.readLatency,            sub->wwlDecoder.writeLatency,
                          sub->wwlDecoder.leakage,                sub->wwlDecoder.area},
      {"wwlSwitchMatrix", SM_R(sub->wwlSwitchMatrix), SM_W(sub->wwlSwitchMatrix),
                          sub->wwlSwitchMatrix.readLatency,       sub->wwlSwitchMatrix.writeLatency,
                          sub->wwlSwitchMatrix.leakage,           sub->wwlSwitchMatrix.area},
      {"wplDecoder",      sub->wplDecoder.readDynamicEnergy,      sub->wplDecoder.writeDynamicEnergy,
                          sub->wplDecoder.readLatency,            sub->wplDecoder.writeLatency,
                          sub->wplDecoder.leakage,                sub->wplDecoder.area},
      {"wplSwitchMatrix", SM_R(sub->wplSwitchMatrix), SM_W(sub->wplSwitchMatrix),
                          sub->wplSwitchMatrix.readLatency,       sub->wplSwitchMatrix.writeLatency,
                          sub->wplSwitchMatrix.leakage,           sub->wplSwitchMatrix.area},
      {"sslDecoder",      sub->sslDecoder.readDynamicEnergy,      sub->sslDecoder.writeDynamicEnergy,
                          sub->sslDecoder.readLatency,            sub->sslDecoder.writeLatency,
                          sub->sslDecoder.leakage,                sub->sslDecoder.area},
      {"sslSwitchMatrix", SM_R(sub->sslSwitchMatrix), SM_W(sub->sslSwitchMatrix),
                          sub->sslSwitchMatrix.readLatency,       sub->sslSwitchMatrix.writeLatency,
                          sub->sslSwitchMatrix.leakage,           sub->sslSwitchMatrix.area},
      {"planeDecoder",    havePlaneDecoder ? planeDecoder.readDynamicEnergy  : 0.0,
                          havePlaneDecoder ? planeDecoder.writeDynamicEnergy : 0.0,
                          havePlaneDecoder ? planeDecoder.readLatency        : 0.0,
                          havePlaneDecoder ? planeDecoder.writeLatency       : 0.0,
                          havePlaneDecoder ? planeDecoder.leakage            : 0.0,
                          havePlaneDecoder ? planeDecoder.area               : 0.0},
    };
    const int NP = (int)(sizeof(P)/sizeof(P[0]));

    if (!isFeRAM) {
        sub->currentSenseAmp.Initialize(numCol, false, false, param->clkFreq, 1);
        sub->currentSenseAmp.CalculateUnitArea();
        sub->currentSenseAmp.CalculateArea(sub->widthArray);
    }

    sub->currentSenseAmp.CalculateLatency(colRes, 1, 1);
    sub->currentSenseAmp.CalculatePower(colRes, 1);
    const double eSenseAmp = sub->currentSenseAmp.readDynamicEnergy;
    //const double tSenseAmp = sub->currentSenseAmp.readLatency;
    
    const double vMinSense = (argVmin > 0) ? argVmin : param->senseAmpResolution;
    const double vAmpIn    = param->senseAmpVmin;

    const double iOffCell   = param->ioffReadTr
                            * pow(10.0, -param->vthReadTr / param->ssReadTr);
    const double iLeakRows  = (double)(numRow - 1) * iOffCell;
    const double iDiffIdeal = MAX(iRead1 - iRead0 - iLeakRows, 1e-15);

    const int    senseMode  = (param->senseMode == 1) ? 1 : 0;
    const double vT         = 8.617333e-5 * (double)param->temp;              
    const double nSub       = param->ssReadTr / (vT * log(10.0));             
    const double senseRefFrac = (param->senseRefFraction > 0.0 && param->senseRefFraction <= 1.0)
                              ? param->senseRefFraction : 0.5;

    const double rRslPath   = twoT ? resCol : 0.0;
    const double vov1       = vIntBias + 0.5*dVmod - param->vthReadTr;        
    const double vov0       = vIntBias - 0.5*dVmod - param->vthReadTr;        
    const double iRead1Eff  = (senseMode == 0)
                            ? DegenerateCurrent(vov1, rRslPath, param->ioffReadTr, param->ssReadTr, param->ionSatReadTr)
                            : iRead1;
    const double iRead0Eff  = (senseMode == 0)
                            ? DegenerateCurrent(vov0, rRslPath, param->ioffReadTr, param->ssReadTr, param->ionSatReadTr)
                            : iRead0;
    const double vSourceIR  = iRead1Eff * rRslPath;
    const double degenLoss  = (iRead1 > 0.0) ? 1.0 - iRead1Eff / iRead1 : 0.0;
    const double iDiffEff   = MAX(iRead1Eff - iRead0Eff - iLeakRows, 1e-15);
    const double iDiff      = (senseMode == 0) ? iDiffEff : iDiffIdeal;      
    const double iSig       = MAX(senseRefFrac * iDiff, 1e-15);              

    const double iAvgCell   = p1 * iRead1 + (1.0 - p1) * iRead0;
    const double resRBL     = resRow * ((param->rblResScale > 0.0) ? param->rblResScale : 1.0);
    const double vDropRBL   = twoT ? (double)numCol * iAvgCell * resRBL * 0.5 : 0.0;
    const double vDsat1     = MAX(vov1, 4.0 * vT);                            
    const double rblDropRatio = (V_RBL_BIAS > vDsat1) ? vDropRBL / (V_RBL_BIAS - vDsat1) : 1e9;

    const double capColSense = twoT ? capRSL : capSSL;
    const double gClamp     = (param->iClampBias > 0.0 ? param->iClampBias : param->ionSatReadTr) / (nSub * vT);
    double tColSettle, tSenseDev, capSenseNode;
    if (senseMode == 0) {
        const double tauCol = capColSense / gClamp + 0.5 * resCol * capColSense;
        tColSettle   = 3.0 * tauCol;
        capSenseNode = MAX(param->capSenseInt, 0.5e-15);
        tSenseDev    = capSenseNode * vAmpIn / iSig;
    } else {
        tColSettle   = 0.5 * resCol * capColSense;
        capSenseNode = capColSense + MAX(param->capSenseAmpIn, 0.0);
        tSenseDev    = capSenseNode * vAmpIn / iSig;
    }

    double tauLatch = 0.0, tSenseRegen = 0.0;
    {
        const double wN = MIN_NMOS_SIZE * tech.featureSize;
        const double wP = tech.pnSizeRatio * MIN_NMOS_SIZE * tech.featureSize;
        double capLatchIn = 0.0, capLatchOut = 0.0;
        CalculateGateCapacitance(INV, 1, wN, wP, tech.featureSize * MAX_TRANSISTOR_HEIGHT,
                                 tech, &capLatchIn, &capLatchOut);
        const double gmLatch = CalculateTransconductance(wN, NMOS, tech)
                             + CalculateTransconductance(wP, PMOS, tech);
        tauLatch = (gmLatch > 0.0)
                 ? (capLatchIn + capLatchOut + MAX(param->capSenseAmpOut, 0.0)) / gmLatch : 0.0;
        const double vTarget = 0.5 * tech.vdd;                 
        tSenseRegen = (tauLatch > 0.0 && vTarget > vAmpIn)
                    ? tauLatch * log(vTarget / vAmpIn) : 0.0;
        if (param->senseAmpRegenTime > 0.0) tSenseRegen = param->senseAmpRegenTime;
        else if (!(tSenseRegen > 5e-12 && tSenseRegen < 2e-9)) {
            cout << "  *** latch regeneration derived as " << tSenseRegen
                 << " s (tau = " << tauLatch << " s) - outside 5 ps..2 ns, using 0.1 ns.\n"
                    "      Set param->senseAmpRegenTime to pin it.\n";
            tSenseRegen = 0.1e-9;
        }
    }
    const double tSenseInt  = tSenseDev;                          
    const double tSenseAmp  = tSenseDev + tSenseRegen;
    const double tSenseNeuroSim = sub->currentSenseAmp.readLatency; 

    const double aSenseAmp = sub->currentSenseAmp.area;
    const double lkSenseAmp= sub->currentSenseAmp.leakage;

    const double LATCH_CLK = 10e9;      // 100 ps clk to Q + setup 
    DFF outLatch(inputParameter, tech, cell);
    outLatch.Initialize(ioWidth, LATCH_CLK);
    outLatch.CalculateArea(0, 0, NONE);
    outLatch.CalculateLatency(1e20, 1);
    outLatch.CalculatePower(1, ioWidth, false);
    const double eLatch = outLatch.readDynamicEnergy;
    const double tLatch = outLatch.readLatency;
    const double aLatch = outLatch.area;
    const double lkLatch= outLatch.leakage;

     // 8.  TIER SUMS
    double eDecRd = 0, eDecWr = 0, aDec = 0, lkDec = 0;
    for (int i = 0; i < NP; i++) {
        eDecRd += P[i].eR; eDecWr += P[i].eW; aDec += P[i].ar; lkDec += P[i].lk;
    }

    const double tColDiffusion = tColSettle;
    const double eReadDCPage = V_RBL_BIAS * iDCPage * (tSenseAmp + tColDiffusion);

          double eReadCore  = eCellReadPage + eReadDCPage + eReadLines + eRestorePage;
          double eWriteCore = eWriteCells + eWriteInhibit + eWriteLines;

          double eReadMacro  = eReadCore  + ePrecharge + eDecRd + eSenseAmp + eLatch;
          double eWriteMacro = eWriteCore + eDecWr;

    double tSelect = 0;
    for (int i = 0; i < NP; i++) tSelect = MAX(tSelect, P[i].tR);
    double tSelectW = 0;
    for (int i = 0; i < NP; i++) tSelectW = MAX(tSelectW, P[i].tW);

    const double kLine     = (param->kLineRc > 0.0) ? param->kLineRc : 1.0;
    const double tPlateRC  = kLine * resRow * capWBL;            

    const double tTransfer = 0.0;                                

    double rPlateEff = rPlateDrv;
    {
        const double tSM = sub->wblSwitchMatrix.readLatency;
        if (tSM > 0.0 && capWBL > 0.0) {
            const double r = (tSM / 0.6931 - resRow * capWBL * 0.5) / capWBL;
            if (r > 1.0 && r < 1.0e7) rPlateEff = r;
        }
    }
    const double rPlateDrive = (param->resPlateDriver > 0.0) ? param->resPlateDriver : rPlateDrv;
    const double tPlateDrive = 2.3 * rPlateDrive * capWBL;       
    const double qPageRead   = (double)numCol * p1 * alphaRd * Qsw;
    const double tPlateSupply= (Vr > vIntBias) ? 2.2 * rPlateDrive * (qPageRead / (Vr - vIntBias)) : 0.0;
    const double tSwitch     = MAX(param->chargeDelay, tPlateSupply);
    const double tCharge     = tPlateDrive + tPlateRC + tSwitch + tTransfer;
    
    const double tSettleRow = kLine * (twoT ? resRBL * capRBL : resRow * capWBL);
    const double tSettleCol = tColDiffusion;
    const double tSettle    = tSettleRow + tSettleCol;

    const double tPre     = 0.0;                    
    const double tRestore = restoreDuty * T_WRITE;  

          double tReadCore  = tCharge + tSettle + tRestore;

          double tReadMacro = tPre + tSelect + tCharge + tSettle
                            + tSenseAmp + tLatch + tRestore;
	const double tLineWrite = 2.3*rPlateDrive*capWBL + kLine*resRow*MAX(capWPL, capWBL);
	  
	const double tPillarWrite = 2.2 * rPlateDrv * capSense;
	      double tWriteCore = nWritePulses * (T_WRITE + tLineWrite + tPillarWrite);
          double tWriteMacro = tSelectW + tWriteCore;

    const double CELL_ENDURANCE   = 1e12;
    const double switchesPerWrite = flipFrac;
    const double writesToFailure  = (flipFrac > 0.0)
                                  ? CELL_ENDURANCE / flipFrac : 1e30;

    const double eWriteTrue = eWriteMacro + (owesARead ? eReadMacro : 0.0);
    const double tWriteTrue = tWriteMacro + (owesARead ? tReadMacro : 0.0);

    const double lkArray = isFeRAM
                         ? (double)numRow * (double)numCol * iOffCell * V_RBL_BIAS
                         : 0.0;
          double lkMacro = lkDec + lkSenseAmp + lkLatch + lkArray;

    const double aArray = sub->areaArray;      
    const double aStair = sub->areaStaircase;  
          double aCore  = aArray + aStair;
          double aMacro = aCore + aDec + aSenseAmp + aLatch;

    double eAddrRd = 0.0, eAddrWr = 0.0, aAddr = 0.0, lkAddr = 0.0;

    if (!isFeRAM) {
        const double aRW = 1.0 / (double)numRow;   

        eReadCore = sub->readDynamicEnergyArray;

        if (cellType == 1) {                      
            eAddrRd = sub->wlDecoder.readDynamicEnergy
                    + sub->precharger.readDynamicEnergy;
            eAddrWr = sub->wlDecoder.writeDynamicEnergy
                    + sub->sramWriteDriver.writeDynamicEnergy;
            aAddr   = sub->wlDecoder.area + sub->precharger.area
                    + sub->sramWriteDriver.area;
            lkAddr  = sub->wlDecoder.leakage + sub->precharger.leakage
                    + sub->sramWriteDriver.leakage;
        } else {                                   
            eAddrRd = sub->wlDecoder.readDynamicEnergy
                    + sub->wlDecoderDriver.readDynamicEnergy
                    + sub->wlNewDecoderDriver.readDynamicEnergy
                    + sub->slSwitchMatrix.readDynamicEnergy;
            eAddrWr = sub->wlDecoder.writeDynamicEnergy
                    + sub->wlDecoderDriver.writeDynamicEnergy
                    + sub->wlNewDecoderDriver.writeDynamicEnergy
                    + sub->slSwitchMatrix.writeDynamicEnergy;
            aAddr   = sub->wlDecoder.area + sub->wlDecoderDriver.area
                    + sub->wlNewDecoderDriver.area + sub->slSwitchMatrix.area;
            lkAddr  = sub->wlDecoder.leakage + sub->wlDecoderDriver.leakage
                    + sub->wlNewDecoderDriver.leakage + sub->slSwitchMatrix.leakage;
        }

        if (cellType == 1) {
            eWriteCore = cell.capSRAMCell * tech.vdd * tech.vdd
                         * (double)numCol * (double)numRow * aRW          
                       + sub->capRow1 * tech.vdd * tech.vdd * (double)numRow * aRW  
                       + sub->capCol  * tech.vdd * tech.vdd * (double)numCol;       
            const double rAccessSRAM = CalculateOnResistance(
                    cell.widthSRAMCellNMOS * tech.featureSize, NMOS,
                    inputParameter.temperature, tech);
            tWriteCore = 2.3 * rAccessSRAM * cell.capSRAMCell;

        } else {
            const double tPulseRRAM = T_WRITE_RRAM;
            const double Vw2   = cell.writeVoltage * cell.writeVoltage;
            const double eCell = 0.5 * Vw2 * tPulseRRAM
                               * (1.0 / cell.resMemCellOn + 1.0 / cell.resMemCellOff);
            eWriteCore = (double)numCol * eCell
                       + sub->capCol * Vw2 * (double)numCol * aRW * (double)numRow
                       + sub->capRow2 * cell.accessVoltage * cell.accessVoltage
                         * (double)numRow * aRW;
            tWriteCore = tPulseRRAM;               
        }

        eReadMacro  = eReadCore  + eAddrRd + eSenseAmp + eLatch;
        eWriteMacro = eWriteCore + eAddrWr;

        tReadCore   = sub->readLatency;
        tReadMacro  = tReadCore + tSenseAmp + tLatch;
        tWriteMacro = tWriteCore + tSelectW
                    + ((cellType == 1) ? sub->sramWriteDriver.writeLatency : 0.0);

        aCore   = sub->areaArray;
        aMacro  = aCore + aAddr + aSenseAmp + aLatch;
        lkMacro = lkAddr + lkSenseAmp + lkLatch;

        if (!(eWriteCore > 0.0) || !(tWriteCore > 0.0) ||
            !(eReadCore  > 0.0) || !(tReadCore  > 0.0)) {
            cout << "\n*** " << cn << ": a tier came back as ZERO"
                 << "  (eReadCore=" << eReadCore << " eWriteCore=" << eWriteCore
                 << " tReadCore=" << tReadCore  << " tWriteCore=" << tWriteCore << ")\n"
                 << "    Do NOT plot this row.  Either the cell branch in "
                 << "SubArray.cpp did not execute\n"
                 << "    (check conventionalSequential) or a parasitic is "
                 << "uninitialised for this cell type.\n";
        }
    }


     // 9.  SENSE-MARGIN DESIGN RULE CHECK
    const double dVsense = dVmod;
    const double vLevel0  = CAP_FE_LIN * Vr / capSense;
    const double vRef     = vLevel0 + 0.5 * dVsense;
    const bool   sensible = isFeRAM ? (dVsense > vMinSense) : true;

     // 10.  REPORT
    const double bits = numCol;
    cout << "\n==============================================================================\n"
         << "  " << cn << "   read-out: " << rn << "   n = " << nPlanes
         << "   " << param->technode << " nm\n"
         << "  READ / WRITE OPERATOR MODE -- single row, full memory periphery,\n"
         << "  CIM periphery excluded\n"
         << "==============================================================================\n";
    H("CONFIGURATION");
    Rf("physical rows (pillars)", (double)numRow, "", 0);
    Rf("columns = page width", (double)numCol, "bits", 0);
    Rf("stacked planes n", (double)nPlanes, "", 0);
    Rf("I/O width", (double)ioWidth, "bits", 0);
    Rf("capacity of this subarray", (double)numRow*numCol*nPlanes/8.0/1024.0, "KiB", 2);
    Rf("read voltage", cell.readVoltage, "V", 3);
    Rf("write voltage", cell.writeVoltage, "V", 3);
    cout << "    data pattern:   p1 (fraction of '1' in the page read) = " << fixed
         << setprecision(3) << p1 << '\n'
         << "                    q1 (fraction of '1' in the word written) = " << q1 << '\n'
         << "    write scheme:   " << sn << '\n'
         << "    flip fraction:  " << fixed << setprecision(4) << flipFrac
         << "   (cells that actually SWITCH, per written page)\n"
         << "    drive fraction: " << fixed << setprecision(4) << driveFrac
         << "   (field applications per cell -- the C*V^2 population)\n"
         << "    plate phases:   " << (int)nWritePulses
         << (owesARead ? "   + a preceding read\n" : "\n");
    if (fabs(cell.readVoltage - cell.writeVoltage) < 1e-12 && cell.mem_rdo != Type::dro) {
        cout << "\n  *** WARNING: readVoltage (" << cell.readVoltage
             << " V) == writeVoltage, but mem_rdo is not 'dro'.\n"
             << "      Param.cpp has a single readVoltage that ignores the read-out mode,\n"
             << "      so this run applies a FULL switching field and is physically a\n"
             << "      destructive read. alpha_rd and the 1/N restore duty below are then\n"
             << "      unjustified. Re-run with a reduced read voltage as argv[11],\n"
             << "      e.g.  ./rw_main 5 2 " << nPlanes << " " << numRow << " " << numCol
             << " 0.5 0.5 2 0 0.025 1.2\n";
    }
    if (cell.mem_rdo == Type::qndro)
        cout << "    qndro:          switches " << param->qndroSwitchFraction*100
             << "% of Pr per read; restore every " << param->qndroRefreshInterval
             << " reads\n";

    H("SENSE MARGIN  (design rule check -- read this before any energy number)");
    cout << "    sensed node: " << senseNodeName << '\n';
    R("Qsw = 2Pr * A_FE", Qsw, "C");
    R("C_sense", capSense, "F");
    R("signal split  dV = Qsw/(C+C_FE)", dVsense, "V");
    Rf("  in millivolts", dVsense * 1e3, "mV", 3);
    R("'0' level (common mode)", vLevel0, "V");
    R("reference level must sit at", vRef, "V");
    R("half margin   dV/2", 0.5*dVsense, "V");
    R("amplifier resolution", vMinSense, "V");
    cout << "    VERDICT: " << (sensible ? "PASS -- '1' and '0' are distinguishable"
                                         : "*** FAIL -- NOT SENSIBLE, energy numbers below are meaningless ***")
         << '\n';
    if (twoT)
        Note("2T-nC is PLANE limited: all n capacitors share the pillar node, so "
             "dV falls as 1/n. Columns do not hurt it.");
    else
        Note("1T-nC is COLUMN limited: the charge lands on the bit line, so "
             "dV falls as 1/numCol. Planes do not hurt it.");
    if (!sensible)
        Note("Fix by: fewer planes per pillar (2T) or fewer columns per bit line "
             "(1T), a larger capacitor area, a segmented/hierarchical line, or a "
             "charge-transfer / offset-cancelled amplifier.");

    H("READ -- energy per access (one page)");
    R("cell array, switching term", numCol * p1 * alphaRd * Qsw * Vr, "J");
    R("cell array, linear term", numCol * eCellRead0, "J");
    if (twoT) R("read-transistor DC path", eReadDCPage, "J");
    R("array lines", eReadLines, "J");
    R("restore (write-back)", eRestorePage, "J");
    R("  = CORE total", eReadCore, "J");
    R("precharge / equalise", ePrecharge, "J");
    R("address decoders + drivers", eDecRd, "J");
    R("sense amplifiers", eSenseAmp, "J");
    R("output latch", eLatch, "J");
    R("  = MACRO total", eReadMacro, "J");
    Rf("  MACRO per bit", eReadMacro / bits * 1e15, "fJ/bit");
    Rf("  periphery share of macro read",
       100.0 * (eReadMacro - eReadCore) / MAX(eReadMacro, 1e-30), "%");

    H("READ -- the 1 versus 0 decomposition (this is the answer to your question)");
    R("energy to read ONE cell storing '1'",
      eCellRead1 + (twoT ? V_RBL_BIAS*iRead1*(tSenseAmp + tColDiffusion) : 0.0), "J");
    R("energy to read ONE cell storing '0'",
      eCellRead0 + (twoT ? V_RBL_BIAS*iRead0*(tSenseAmp + tColDiffusion) : 0.0), "J");
    {
        double e1 = eCellRead1 + (twoT ? V_RBL_BIAS*iRead1*(tSenseAmp + tColDiffusion) : 0.0);
        double e0 = eCellRead0 + (twoT ? V_RBL_BIAS*iRead0*(tSenseAmp + tColDiffusion) : 0.0);
        Rf("ratio  E(1)/E(0)", (e0 > 0 ? e1/e0 : 0.0), "x");
    }
    cout << "\n    page read energy as a function of p1 (everything else fixed):\n"
         << "      p1      cell+DC (J)        restore (J)        MACRO (J)     fJ/bit\n";
    for (double pp = 0.0; pp <= 1.0001; pp += 0.25) {
        double eCell = numCol * (pp*eCellRead1 + (1-pp)*eCellRead0);
        double eDC   = twoT ? numCol * (pp*V_RBL_BIAS*iRead1*(tSenseAmp + tColDiffusion)
                                      + (1-pp)*V_RBL_BIAS*iRead0*(tSenseAmp + tColDiffusion)) : 0.0;
        double eRes  = restoreDuty * (pp*numCol*eFlipCell + eInhibit);
        double eMac  = eCell + eDC + eReadLines + eRes + ePrecharge
                     + eDecRd + eSenseAmp + eLatch;
        cout << "    " << fixed << setprecision(2) << setw(6) << pp
             << scientific << setprecision(4)
             << setw(18) << (eCell+eDC) << setw(19) << eRes << setw(18) << eMac
             << fixed << setprecision(2) << setw(11) << eMac/bits*1e15 << '\n';
    }

    H("WRITE -- energy per access (one page)");
    R("cells (flipFrac*Qsw*V + driveFrac*C*V^2)", eWriteCells, "J");
    R("V/3 inhibition (half-select)", eWriteInhibit, "J");
    R("array lines", eWriteLines, "J");
    R("  = CORE total", eWriteCore, "J");
    R("address decoders + write drivers", eDecWr, "J");
    R("  = MACRO total", eWriteMacro, "J");
    Rf("  MACRO per bit", eWriteMacro / bits * 1e15, "fJ/bit");
    Rf("  inhibition share of macro write",
       100.0 * eWriteInhibit / MAX(eWriteMacro, 1e-30), "%");
    Rf("  switches per cell per write", switchesPerWrite, "");
    Rf("  writes to failure (1e12 cycle cell)", writesToFailure, "");
    if (owesARead) {
        R("  + the read this scheme owes", eReadMacro, "J");
        R("  = TRUE cost per write operation", eWriteTrue, "J");
    }
    cout << "\n    page write energy as a function of the data written:\n"
         << "      case                     flip   drive    cells (J)      MACRO (J)   fJ/bit\n";
    {
        struct C { const char *nm; double pp, qq; } cs[] = {
            {"all-0 over all-0",   0.0, 0.0},
            {"all-0 -> all-1",     0.0, 1.0},
            {"all-1 -> all-0",     1.0, 0.0},
            {"all-1 over all-1",   1.0, 1.0},
            {"random over random", 0.5, 0.5},
        };
        for (int i = 0; i < 5; i++) {
            double hm = cs[i].pp*(1-cs[i].qq) + cs[i].qq*(1-cs[i].pp);
            double ff = (scheme == 2) ? (cs[i].pp + cs[i].qq) : hm;
            double df = (scheme == 2) ? 2.0 : (scheme == 3) ? 1.0 : hm;
            double ec = numCol * (ff*eSwitchCell + df*eLinearCell);
            double em = ec + eWriteInhibit + eWriteLines + eDecWr;
            cout << "    " << left << setw(24) << cs[i].nm << right
                 << fixed << setprecision(3) << setw(7) << ff << setw(8) << df
                 << scientific << setprecision(4) << setw(15) << ec << setw(15) << em
                 << fixed << setprecision(2) << setw(10) << em/bits*1e15 << '\n';
        }
        cout << "\n    (flip = cells that switch -> endurance; drive = field applications\n"
             << "     per cell -> the C*V^2 term.  erase-program is the only scheme whose\n"
             << "     flip count does NOT fall as the data changes less.)\n";
    }

    H("LATENCY -- worst case, therefore DATA INDEPENDENT");
    R("precharge", tPre, "s");
    R("address decode + line drive", tSelect, "s");
    R("charge transfer", tCharge, "s");
    R("  of which FE switching floor", param->chargeDelay, "s");
    R("  of which plate-driver supply bound", tPlateSupply, "s");
    R("  of which plate strip RC (2.2RC)", tPlateRC, "s");
    R("  (removed) R_plate*n*C_FE transfer", 2.2 * rPlateDrv * capSense, "s");
    R("line settle (RBL far end + column)", tSettle, "s");
    R("  column leg, worst-case row", tSettleCol, "s");
    if (senseMode == 0) {
        R("    = 3*(C_col/g_clamp + 0.5*R*C), C_col", capColSense, "F");
        R("      g_clamp", gClamp, "S");
    } else {
        R("    = 0.5*R_col*C_col (far-row diffusion)", 0.5 * resCol * capColSense, "s");
    }
    R("sense amplifier (= develop + regenerate)", tSenseAmp, "s");
    R("  develop 60 mV: C_int*Vmin/(f*dI)", tSenseDev, "s");
    R(senseMode == 0 ? "    C_int (amplifier internal node)" : "    C_int (column + amp input)", capSenseNode, "F");
    R("    usable dI = f_ref*(I1-I0-Ileak)", iSig, "A");
    R("    I(read 1) at the far row (IR-degenerated)", iRead1Eff, "A");
    R("    source IR drop of the far row", vSourceIR, "V");
    Rf("    far-row '1' current loss", 100.0 * degenLoss, "%", 2);
    R("  latch regeneration tau*ln(Vdd/2/Vmin)", tSenseRegen, "s");
    R("    tau_latch = C/gm", tauLatch, "s");
    R("  NeuroSim CurrentSenseAmp (not used)", tSenseNeuroSim, "s");
    R("output latch (clk-to-Q)", tLatch, "s");
    R("restore (amortised)", tRestore, "s");
    R("  = READ, CORE", tReadCore, "s");
    R("  = READ, MACRO (tAA)", tReadMacro, "s");
    R("  = WRITE, CORE", tWriteCore, "s");
    R("  = WRITE, MACRO", tWriteMacro, "s");
    Rf("  read bandwidth", bits / tReadMacro / 1e9, "Gb/s");

    H("AREA AND LEAKAGE");
    R("cell array", sub->areaArray, "m^2");
    R("staircase contacts", sub->areaStaircase, "m^2");
    R("  = CORE", aCore, "m^2");
    R("decoders + switch matrices", aDec, "m^2");
    R("sense amplifiers", aSenseAmp, "m^2");
    R("output latch", aLatch, "m^2");
    R("  = MACRO", aMacro, "m^2");
    Rf("  array efficiency (core/macro)", 100.0*aCore/MAX(aMacro,1e-30), "%");
    Rf("  density", (double)numRow*numCol*nPlanes/(aMacro*1e12), "Mbit/mm^2", 2);
    R("macro leakage", lkMacro, "W");
    R("resRow / resCol", resRow, "ohm");
    R("  resCol", resCol, "ohm");

    H("PER-BLOCK DETAIL  (included in the MACRO tier)");
    cout << "  " << left << setw(20) << "block"
         << right << setw(15) << "read E (J)" << setw(15) << "write E (J)"
         << setw(14) << "read t (s)" << setw(14) << "area (m^2)" << '\n';
    for (int i = 0; i < NP; i++)
        cout << "  " << left << setw(20) << P[i].name << right << scientific
             << setprecision(4) << setw(15) << P[i].eR << setw(15) << P[i].eW
             << setw(14) << P[i].tR << setw(14) << P[i].ar << '\n';
    cout << "  " << left << setw(20) << "precharge" << right << scientific
         << setprecision(4) << setw(15) << ePrecharge << setw(15) << 0.0
         << setw(14) << tPre << setw(14) << 0.0 << '\n';
    cout << "  " << left << setw(20) << "currentSenseAmp" << right << scientific
         << setprecision(4) << setw(15) << eSenseAmp << setw(15) << 0.0
         << setw(14) << tSenseAmp << setw(14) << aSenseAmp << '\n';
    cout << "  " << left << setw(20) << "output latch (dff)" << right << scientific
         << setprecision(4) << setw(15) << eLatch << setw(15) << 0.0
         << setw(14) << tLatch << setw(14) << aLatch << '\n';
    cout << "  " << left << setw(20) << "[cell array]" << right << scientific
         << setprecision(4) << setw(15) << (eReadCore) << setw(15) << eWriteCore
         << setw(14) << tReadCore << setw(14) << aCore << '\n';

    H("EXCLUDED  (CIM periphery -- printed so the exclusion is auditable)");
    cout << "  " << left << setw(24) << "sarADC"
         << right << scientific << setprecision(4) << setw(15)
         << sub->sarADC.readDynamicEnergy << "  J   (not in MACRO)\n";
    cout << "  " << left << setw(24) << "multilevelSenseAmp"
         << right << setw(15) << sub->multilevelSenseAmp.readDynamicEnergy << "  J\n";
    cout << "  " << left << setw(24) << "multilevelSAEncoder"
         << right << setw(15) << sub->multilevelSAEncoder.readDynamicEnergy << "  J\n";
    cout << "  " << left << setw(24) << "mux + muxDecoder"
         << right << setw(15) << (sub->mux.readDynamicEnergy + sub->muxDecoder.readDynamicEnergy) << "  J\n";
    cout << "  " << left << setw(24) << "shiftAddWeight"
         << right << setw(15) << sub->shiftAddWeight.readDynamicEnergy << "  J\n";
    cout << "  " << left << setw(24) << "shiftAddInput"
         << right << setw(15) << sub->shiftAddInput.readDynamicEnergy << "  J\n";
    cout << "  " << left << setw(24) << "adder (accumulation)"
         << right << setw(15) << sub->adder.readDynamicEnergy << "  J\n";
    cout << "  " << left << setw(24) << "beta (CIM validation)"
         << right << fixed << setprecision(2) << setw(15) << param->beta << "      not applied\n";

    H("CORRECTIONS APPLIED (FIX C and FIX G)");
    {
        double dffTotal = sub->rblSwitchMatrix.dff.readDynamicEnergy
                        + sub->rslSwitchMatrix.dff.readDynamicEnergy
                        + sub->wblSwitchMatrix.dff.readDynamicEnergy
                        + sub->wwlSwitchMatrix.dff.readDynamicEnergy
                        + sub->wplSwitchMatrix.dff.readDynamicEnergy
                        + sub->sslSwitchMatrix.dff.readDynamicEnergy;
        R("scan-chain DFF energy removed", dffTotal, "J");
        R("  of which wblSwitchMatrix", sub->wblSwitchMatrix.dff.readDynamicEnergy, "J");
        Note("A SwitchMatrix holds one DFF per output so a CIM activation vector can");
        Note("be scanned in. A memory macro selects the line from the decoder and has");
        Note("no such register, so this is hardware the design does not build.");
    }
    R("wblSwitchMatrix, raw from SubArray", sub->wblSwitchMatrix.readDynamicEnergy, "J");
    R("wblSwitchMatrix, scan chain removed", SM_R(sub->wblSwitchMatrix), "J");
    R("wblSwitchMatrix, per ONE activation",
      SM_R(sub->wblSwitchMatrix) * WBL_ACT_FIX, "J");
    Rf("correction factor", 1.0/WBL_ACT_FIX, "x  (= bitsPerCell)", 1);
    R("rblSwitchMatrix, scan chain removed", SM_R(sub->rblSwitchMatrix), "J");
    Note("These last two must be IDENTICAL. A SwitchMatrix's read energy is its own");
    Note("TG drain and gate caps only -- the LINE capacitance is not inside it (it is");
    Note("counted once in 'array lines' above), so one activation of a WBL driver and");
    Note("one activation of an RBL driver cost the same. If they differ, the");
    Note("correction is the wrong shape: set WBL_ACT_FIX to 1.0 and investigate");
    Note("SwitchMatrix::CalculatePower directly.");

    H("CROSS-CHECK AGAINST SubArray.cpp");
    {
        double saRead  = sub->readDynamicEnergyArray;
        double saWrite = sub->writeDynamicEnergyArray;
        R("SubArray readDynamicEnergyArray", saRead, "J");
        R("SubArray writeDynamicEnergyArray", saWrite, "J");
        R("this driver, core read", eReadCore, "J");
        R("this driver, core write", eWriteCore, "J");
        Note("These differ on purpose. SubArray.cpp hard-codes a 50/50 pattern and");
        Note("alpha_Fe = 0.5 (a differential write on random data); it also folds");
        Note("write-back PERIPHERY into readDynamicEnergyArray. Set p1 = q1 = 0.5 and");
        Note("scheme = 1 to compare the cell terms like for like.");
    }

    H("STRUCTURAL CHECKS");
    cout << "  numRow * activityRowRead        = " << fixed << setprecision(6)
         << numRow * sub->activityRowRead  << "   (must be 1.000000)\n";
    cout << "  numRow * activityRowWrite       = "
         << numRow * sub->activityRowWrite << "   (must be 1.000000)\n";
    cout << "  core read  < macro read         : " << ((eReadCore  < eReadMacro)  ? "ok" : "FAIL") << '\n';
    cout << "  core write < macro write        : " << ((eWriteCore < eWriteMacro) ? "ok" : "FAIL") << '\n';
    cout << "  ndro restore == 0               : "
         << ((cell.mem_rdo != Type::ndro || eRestorePage == 0.0) ? "ok" : "FAIL") << '\n';
    cout << "  all-0 read cheaper than all-1   : "
         << ((eCellRead0 <= eCellRead1) ? "ok" : "FAIL") << '\n';
    cout << "  every energy > 0                : "
         << ((eReadMacro > 0 && eWriteMacro > 0) ? "ok" : "FAIL -- a term is still commented out") << '\n';
    cout << "  sense margin                    : " << (sensible ? "ok" : "FAIL") << '\n';
    if (cellType == 2)
        cout << "  RRAM R_on / R_off               : " << cell.resistanceOn << " / "
             << cell.resistanceOff << " ohm   (FIX A: re-derived, NOT Param's "
             << "memcelltype==5 value of 1e20)\n";
    cout << "  plate-driver resistance         : " << rPlateDrv << " ohm  "
         << (rPlateIsReal ? "(wblSwitchMatrix.resTg -- used for tPillarWrite only)"
                          : "(FALLBACK 5 kOhm placeholder -- resTg was not usable)") << '\n';
    cout << "  plate driver NeuroSim timed     : " << rPlateEff
         << " ohm  (back-solved from wblSwitchMatrix.readLatency; sets tPlateSupply)\n";
    cout << "  plate supply bound vs FE floor  : " << tPlateSupply << " s vs "
         << param->chargeDelay << " s  " << (tPlateSupply > param->chargeDelay ? "DRIVER-LIMITED" : "FE-limited") << '\n';
    cout << "  sense: develop / regen / total  : " << tSenseDev << " / " << tSenseRegen
         << " / " << tSenseAmp << " s   (NeuroSim block said " << tSenseNeuroSim << " s)\n";
    cout << "  read lines: (n-1) unselected strips " << eReadUnselStrips << " J of " << eReadLines << " J\n";
    cout << "  sense mode                      : " << (senseMode == 0 ? "0 = clamped column (current-mode)" : "1 = column integrates (voltage-mode)") << '\n';
    cout << "  far-row source IR degeneration  : " << vSourceIR << " V, '1' current -" << 100.0*degenLoss << " %\n";
    cout << "  row-line (RBL) IR drop, page    : " << vDropRBL << " V  ratio to (V_RBL - V_dsat) = " << rblDropRatio
         << (rblDropRatio > 1.0 ? "  *** PAGE CURRENT COLLAPSES V_DS AT THE FAR COLUMNS -- widen/strap RBL or read a narrower page ***" : "  ok") << '\n';
    cout << "  signal clipped to rail          : " << (dVclipped ? "YES -- dV limited by (Vr - vIntBias)" : "no") << '\n';

    cerr << cn << ',' << rn << ',' << nPlanes << ',' << numRow << ',' << numCol << ','
         << p1 << ',' << q1 << ','
         << (scheme == 2 ? "erase-program" : scheme == 3 ? "vanilla" : "differential") << ','
         << flipFrac << ','
         << eReadCore  << ',' << eReadMacro  << ',' << eReadMacro/bits*1e15 << ','
         << eWriteCore << ',' << eWriteMacro << ',' << eWriteMacro/bits*1e15 << ','
         << tReadCore  << ',' << tReadMacro  << ',' << tWriteCore << ',' << tWriteMacro << ','
         << aCore*1e12 << ',' << aMacro*1e12 << ',' << lkMacro << ','
         << aArray*1e12 << ',' << aStair*1e12 << ','
         << dVsense << ','
         << (!sensible ? "FAIL" : dVclipped ? "PASS (clipped)" : "PASS") << ','
         << (isFeRAM ? eDecRd : eAddrRd) << ','      
         << (isFeRAM ? eDecWr : eAddrWr) << ','
         << eSenseAmp << ',' << ePrecharge << ',' << eLatch;

    cerr << ',' << (isFeRAM ? CAP_FE_LIN : 0.0)
         << ',' << (isFeRAM ? capSense   : 0.0)
         << ',' << (isFeRAM ? Qsw        : 0.0)
         << ',' << (isFeRAM ? vIntBias   : 0.0)
         << ',' << (isFeRAM ? vUnsel     : 0.0)
         << ',' << (isFeRAM ? iRead1     : 0.0)
         << ',' << (isFeRAM ? iRead0     : 0.0)
         << ',' << (isFeRAM ? T_WRITE    : tWriteCore)
         << ',' << (isFeRAM ? nWritePulses : 1.0);

    cerr << ',' << (isFeRAM ? numCol * p1 * alphaRd * Qsw * Vr : 0.0)
         << ',' << (isFeRAM ? numCol * eCellRead0              : 0.0)
         << ',' << (isFeRAM ? eReadDCPage                      : 0.0)
         << ',' << (isFeRAM ? eReadLines                       : eReadCore)
         << ',' << (isFeRAM ? eRestorePage                     : 0.0);

    cerr << ',' << (isFeRAM ? eWriteCells   : eWriteCore)
         << ',' << (isFeRAM ? eWriteInhibit : 0.0)
         << ',' << (isFeRAM ? eWriteLines   : 0.0)
	 << ',' << driveFrac << ',' << switchesPerWrite << ',' << writesToFailure
         << ',' << eWriteTrue << ',' << tWriteTrue;

    cerr << ',' << nHalfSel
         << ',' << vInh << ',' << eHalfMVcm << ',' << disturbRatio;

    cerr << ',' << tSelect << ',' << tCharge << ',' << tSettle
         << ',' << tSenseAmp << ',' << tLatch << ',' << tRestore << ',' << tSelectW;
    cerr << ',' << tSwitch << ',' << tPlateRC << ',' << tTransfer
         << ',' << tSettleRow << ',' << tSettleCol;
    cerr << ',' << tSenseInt << ',' << iDiff << ',' << iLeakRows
         << ',' << capSenseNode << ',' << readOv;

    for (int i = 0; i < NP; i++)
        cerr << ',' << P[i].eR << ',' << P[i].eW << ',' << P[i].tR
             << ',' << P[i].tW << ',' << P[i].lk << ',' << P[i].ar;
    cerr << ',' << eSenseAmp   << ",0,"  << tSenseAmp << ",0," << lkSenseAmp << ',' << aSenseAmp
         << ',' << eLatch      << ",0,"  << tLatch    << ",0," << lkLatch    << ',' << aLatch
         << ',' << eReadCore   << ','    << eWriteCore << ',' << tReadCore
         << ',' << tWriteCore  << ','    << lkArray    << ',' << aCore;

    cerr << ',' << tSenseDev << ',' << tSenseRegen << ',' << tauLatch
         << ',' << tSenseNeuroSim << ',' << tPlateSupply
         << ',' << rPlateDrv << ',' << rPlateEff
         << ',' << param->capJunctionTr << ',' << capWBLpar << ',' << capWBLwire
         << ',' << eReadUnselStrips << ',' << eWriteUnselStrips << ',' << senseRefFrac;
    cerr << ',' << senseMode << ',' << gClamp << ',' << capColSense
         << ',' << iRead1 << ',' << iRead0 << ',' << rRslPath
         << ',' << vSourceIR << ',' << degenLoss << ',' << vDropRBL << ',' << rblDropRatio
         << ',' << kLine;
    cerr << '\n';

    delete sub;
    return sensible ? 0 : 2;
}
