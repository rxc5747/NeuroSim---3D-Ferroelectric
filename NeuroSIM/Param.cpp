/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
*   
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer. 
*   
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
*   
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cstdio>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "math.h"
#include "Param.h"
#include "typedef.h"
#include "constant.h"

using namespace std;

Param::Param() {
	/***************************************** user defined design options and parameters *****************************************/
	operationmode = 2;     		// 1: conventionalSequential (Use several multi-bit RRAM as one synapse)
								// 2: conventionalParallel (Use several multi-bit RRAM as one synapse)

	memcelltype = 5;        	// 1: cell.memCellType = Type::SRAM
								// 2: cell.memCellType = Type::RRAM
								// 3: cell.memCellType = Type::FeFET
								// 4: cell.memCellType = Type::nvCap
								// 5: cell.memCellType = Type::_2TnC
								// 6: cell.memCellType = Type::_1TnC
								// 7: cell.memCellType = Type::_1T1C
	
	mem_rdo = 3;         	// 1: cell.mem_rdo = Type::ndro
								// 2: cell.mem_rdo = Type::qndro
								// 3: cell.mem_rdo = Type::dro

	
	accesstype = 1;         	// 1: cell.accessType = CMOS_access
								// 2: cell.accessType = BJT_access
								// 3: cell.accessType = diode_access
								// 4: cell.accessType = none_access (Crossbar Array)
	
	transistortype = 1;     	// 1: inputParameter.transistorType = conventional
	
	deviceroadmap = 2;      	// 1: inputParameter.deviceRoadmap = HP
								// 2: inputParameter.deviceRoadmap = LSTP
								
	// Anni update
	globalBusType = false;		// false: X-Y Bus
								// true: H-Tree
								
	globalBufferType = false;    // false: register file
								// true: SRAM
	globalBufferCoreSizeRow = 128;
	globalBufferCoreSizeCol = 128;
	
	tileBufferType = false;      // false: register file
								// true: SRAM
	tileBufferCoreSizeRow = 32;
	tileBufferCoreSizeCol = 32;
	
	peBufferType = false;        // false: register file
								// true: SRAM
	
	chipActivation = true;      // false: activation (reLu/sigmoid) inside Tile
								// true: activation outside Tile
						 		
	reLu = true;                // false: sigmoid
								// true: reLu
								
	novelMapping = true;        // false: conventional mapping (change to false for swin_t)
								// true: novel mapping
								
	SARADC = true;              // false: MLSA
	                            // true: sar ADC
	currentMode = false;         // false: MLSA use VSA
	                            // true: MLSA use CSA
								// use VSA for nvCap
	
	pipeline = false;            // false: layer-by-layer process --> huge leakage energy in HP
								// true: pipeline process
	speedUpDegree = 1;          // 1 = no speed up --> original speed
								// 2 and more : speed up ratio, the higher, the faster
								// A speed-up degree upper bound: when there is no idle period during each layer --> no need to further fold the system clock
								// This idle period is defined by IFM sizes and data flow, the actual process latency of each layer may be different due to extra peripheries
	
	validated = true;			// false: no calibration factors
								// true: validated by silicon data (wiring area in layout, gate switching activity, post-layout performance drop...)
								
	synchronous = true;			// false: asynchronous
								// true: synchronous, clkFreq will be decided by sensing delay
	
	memoryMode           = false;
	coreOnly             = false;
	numRowActivated      = 1;
	qndroRefreshInterval = 1e6;     // MUST MATCH ProcessingUnit.cpp 

	// twoPr = 0.20;        // 20 uC/cm^2 Adnaan et al. JXCDC 2025, Fig. 1(a)

	epsFE         = 30.0;      // HZO relative permittivity
	tFE           = 5e-9;     // ALD 10 nm HZO, vertical 2T-nC demo 
	twoPr         = 0.40;      // 40 uC/cm^2, vertical cliff cell 
	ecFE          = 1.2;       // coercive field MV/cm, HZO typical 1.0-1.5
	nlsTauInf     = 236e-9;    // HARD floor - no field switches faster
	nlsEa         = 2.42;      // MV/cm
	nlsAlpha      = 3.73;

	vthReadTr     = 0.30;      // 22 nm LP NMOS
	ssReadTr      = 0.064;     // 64 mV/dec, Intel 22 nm tri-gate NMOS 
	ionSatReadTr  = 10e-6;     // >10 uA measured 
	ioffReadTr    = 10e-9;
	capGateRatio  = 3.0;       // C_MFM ~ 3x C_gate, QNRO design point 
	senseMarginVolt = 0.05;    // 50 mV guard band either side of Vth

	plateSheetRes = 12.5;

	staircaseStep    = 400e-9;   
	staircasePadPitch = 400e-9;  

	inhibitPillarModel = 0;
	inhibitSupplyDraw  = true;
	
	resPlateDriver = 5.0e3;
	
	dataOnesRead  = 0.5;       // 0.5 reproduces the old 50/50 assumption

	readOverdriveAuto  = true;
	readOverdrive      = 0.1375;

	senseAmpResolution = 0.025;

	senseAmpVmin       = 0.060;

	plateContactSides = 2.0;

	senseMode          = 0;

	capSenseInt        = 2.0e-15;

	iClampBias         = 10e-6;

	capJunctionTr      = 0.05e-15;
	capSenseAmpIn      = 0.3e-15;
	capSenseAmpOut     = 2.0e-15;
	senseRefFraction   = 0.5;
	senseAmpRegenTime  = 0.0;
	readDrivesUnselPlanes = true;
	kLineRc            = 1.0;
	rblResScale        = 1.0;

	rhoPillar    = 1.0e-7;    // W-filled pillar (TiN liner + W core)
	pillarRadius = 40e-9;     // = rString in SubArray.cpp
	planePitch   = 60e-9;     // = tWBL + distWBL in SubArray.cpp
	rhoVia       = 1.0e-7;    // W staircase via
	viaRadius    = 40e-9;

	plateDriverSized  = 0;
	plateDriverFanout = 4.0;

	algoWeightMax = 1;
	algoWeightMin = -1;

	inhibitDivider      = 2.0;   // V/2 self-boosted. V/3 measured to FAIL   
	qndroSwitchFraction = 0.20;  // MEASURE from P-V loop at V_read:
	                             // it is dP(V_read)/2Pr. Placeholder.       
	
	clkFreq = 1e9;                      // Clock frequency
	temp = 300;                         // Temperature (K)

	technode = 22;					    // Technology node (nm)
	
	// 1.4 update: Activation implementation option added
	Activationtype=true; // true: SRAM, False: RRAM

	// 1.4 update: special layout (for GAA nodes (1, 2, nm) only)
	speciallayout=1;

	// 1.4 update
	sizingfactor_MUX=1; // sizing for the final driver of mux in rowdecoder.cpp (important for technology scaling)
	sizingfactor_WLdecoder=1; // sizing for the final driver of WLdecoder in rowdecoder.cpp (important for technology scaling)

	// 1.4 update: switchmatrix parameter tuning
	newswitchmatrixsizeratio=6;
	switchmatrixsizeratio=0.1;

	// 1.4 update: buffer parameters : change it if you have to insert buffers
	buffernumber=0;
	buffersizeratio=0;		

	// 1.4 update: new technology node added
	// recommended buffer/driver sizings for MUX, WLdecoder, switchmatrix are provided for each technology node
	// sizingfactor_WLdecoder needs to be adjusted to reduce SRAM WL delay

	switch (technode){ 
		case 22:	sizingfactor_MUX=110; newswitchmatrixsizeratio=14; switchmatrixsizeratio=0.2; break; 
		case 14:	sizingfactor_MUX=130; switchmatrixsizeratio=0.2;  break;  
		case 10:	sizingfactor_MUX=80;  switchmatrixsizeratio=0.2;  break;  
		case 7:		sizingfactor_MUX=60;  switchmatrixsizeratio=0.2;  break;  
		case 5:		sizingfactor_MUX=50;  switchmatrixsizeratio=0.06; break;  
		case 3:		sizingfactor_MUX=25;  switchmatrixsizeratio=0.03; break; 
		case 2:		sizingfactor_MUX=25;  switchmatrixsizeratio=0.1; buffernumber=3; buffersizeratio=3; break;  
		case 1:		sizingfactor_MUX=30;  switchmatrixsizeratio=0.1; buffernumber=3; buffersizeratio=3; break; 
	} 

	// 1.4 update: new wire width, barrierthickness

	switch (technode){
		case 130: 	Metal0=175; Metal1=175; wireWidth=175; barrierthickness=10.0e-9; featuresize = wireWidth*1e-9; break;  
		case 90: 	Metal0=110; Metal1=110; wireWidth=110; barrierthickness=10.0e-9; featuresize = wireWidth*1e-9; break;  
		case 65:	Metal0=105; Metal1=105; wireWidth=105; barrierthickness=7.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 45:	Metal0=80; Metal1=80;   wireWidth=80;  barrierthickness=5.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 32:	Metal0=56; Metal1=56;   wireWidth=56;  barrierthickness=4.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 22:	Metal0=40; Metal1=40;   wireWidth=40;  barrierthickness=2.5e-9;  featuresize = wireWidth*1e-9; break; 
		case 14:	Metal0=32; Metal1=39;   wireWidth=32;  barrierthickness=2.5e-9;  featuresize = wireWidth*1e-9; break;  
		case 10:	Metal0=22; Metal1=32;   wireWidth=22;  barrierthickness=2.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 7:		Metal0=20; Metal1=28.5; wireWidth=20;  barrierthickness=2.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 5:		Metal0=15; Metal1=17;   wireWidth=15;  barrierthickness=2.0e-9;  featuresize = wireWidth*1e-9; break;  
		case 3:		Metal0=12; Metal1=16;   wireWidth=12;  barrierthickness=1.5e-9;  featuresize = wireWidth*1e-9; break; 
		case 2:		Metal0=10; Metal1=11.5; wireWidth=10;  barrierthickness=0.5e-9;  featuresize = wireWidth*1e-9; break;  
		case 1:		Metal0=8;  Metal1=10;   wireWidth=8;   barrierthickness=0.2e-9;  featuresize = wireWidth*1e-9; break;    
		case -1:	break;	
		default:	exit(-1); puts("Wire width out of range"); 
	}
	

	globalBusDelayTolerance = 0.1;      // to relax bus delay for global H-Tree (chip level: communication among tiles), if tolerance is 0.1, the latency will be relax to (1+0.1)*optimalLatency (trade-off with energy)
	localBusDelayTolerance = 0.1;       // to relax bus delay for global H-Tree (tile level: communication among PEs), if tolerance is 0.1, the latency will be relax to (1+0.1)*optimalLatency (trade-off with energy)
	treeFoldedRatio = 4;                // the H-Tree is assumed to be able to folding in layout (save area)
	maxGlobalBusWidth = 2048;           // the max buswidth allowed on chip level (just a upper_bound, the actual bus width is defined according to the auto floorplan)

	// 1.4 update

	inputtoggle = 0.18; // toggling rate for the interconnect energy estimation, needs to be set based on real traces of the input bits

	// Toggling rate tracking can be done by modifying the ProcessingUnit.cpp

	// VGG8-CIFAR10: 0.18
	// ResNet34-ImageNet: 0.29
	// ResNet18-ImageNet: 0.33

	outputtoggle = 0.5; // output bit toggling has a negligible portion of the interconnect energy. Set it to 50 % for simpliciity and generalizability for all neural network workloads.

	if (memcelltype == 5 || memcelltype == 6) {
		numRowSubArrayPhysical = 512;	    // Actual # of rows in single subArray
        	bitsPerCell = 16;		    // # of capacitors ('n' in 2TnC)
        	numRowSubArray = numRowSubArrayPhysical * bitsPerCell;		

		numColSubArray = 512;		    // # of columns in single subArray

		//twoPr      = 0.20;                             // 20 uC/cm^2 Adnaan et al., JXCDC 2025, Fig. 1(a)
		cellAreaFE = 2.0 * 3.14159 * 45e-9 * 40e-9;    // 2*pi*r_s*t_WBL = 1.0053e-14 m^2

	} else {
		numRowSubArrayPhysical = 512;
		numRowSubArray = 512;               // # of rows in single subArray
        	numColSubArray = 512;               // # of columns in single subArray
	}

	// 3D floorplan: architecture + staircase 
	integrationMode = CBA;       // CNA (0) default | CUA (1) | CBA (2)
	staircaseStepWidth     = 200e-9;    // floor on the step run; via geometry may raise it
	viaDiameter3D          = 80e-9;     // = 2 x 40 nm via radius used in Initialize()
	viaOverlayMargin       = 20e-9;
	viaSpacing3D           = 60e-9;
	metalPitch3D           = 100e-9;
	logicPackingEfficiency = 0.70;
	bondPadPitch           = 500e-9;    // shrink as hybrid-bond process allows

	staircaseStepPitch = 400e-9;            // per-plane contact step (NAND: ~0.3-1 um)
	staircaseDummySteps = 2;                // trim/etch margin steps
	staircaseBothSidesFullContact = true;   // every plane contacted both sides (plate RC/2)
	staircaseEdgeMargin = 0.5e-6;           // slit + guard at strip end
	tavPitch = 200e-9;                      // through-array via pitch (CUA)
	cuaUtilization = 0.75;                  // routable fraction under array
	cuaAreaDerate = 1.2;                    // thermal-budget CMOS penalty (CUA only)
	cbaUtilization = 0.85;                  // logic-die placement utilization
	
	numBondPadPerSubarray = numColSubArray + numRowSubArrayPhysical + bitsPerCell + 64;
	
	// 230920 update

	sync_data_transfer=0;

	// 1.4 update: move the operaion mode up

	conventionalParallel = 0;
	conventionalSequential = 0;
	BNNparallelMode = 0;                
	BNNsequentialMode = 0;              
	XNORsequentialMode = 0;          
	XNORparallelMode = 0;         
	switch(operationmode) {
		case 6:	    XNORparallelMode = 1;               break;     
		case 5:	    XNORsequentialMode = 1;             break;     
		case 4:	    BNNparallelMode = 1;                break;     
		case 3:	    BNNsequentialMode = 1;              break;     
		case 2:	    conventionalParallel = 1;           break;     
		case 1:	    conventionalSequential = 1;         break;     
		default:	printf("operationmode ERROR\n");	exit(-1);
	}
	
	/*** parallel read ***/
	parallelRead = 0;
	if(conventionalParallel || BNNparallelMode || XNORparallelMode) {
		parallelRead = 1;
	} else {
		parallelRead = 0;
	}
	// Anni update
	if (parallelRead) {
		if (memcelltype == 5 || memcelltype == 6) {
			numRowParallel = numRowSubArrayPhysical;        // user defined number: >1 and <=numRowSubArray
		} else {
			numRowParallel = numRowSubArrayPhysical;        // user defined number: >1 and <=numRowSubArray
		}

	} else {
		numRowParallel = 1;
	}

	/*** option to relax subArray layout ***/

	relaxArrayCellHeight = 0;           // relax ArrayCellHeight or not
	relaxArrayCellWidth = 0;            // relax ArrayCellWidth or not

	numColMuxed = 8;                    // How many columns share 1 ADC (for eNVM and FeFET) or parallel SRAM
	
	// 1.4 update: handle the exception for conventionalsequential case
	// 1.4 update 230615
	if ((conventionalSequential == 1 || conventionalSequential == 3 || conventionalSequential == 5 ) && (memcelltype==1))
	{
	numColMuxed=numColPerSynapse;
	}
	
	 levelOutput = 32;                   // # of levels of the multilevelSenseAmp output, should be in 2^N forms; e.g. 32 levels --> 5-bit ADC
         // levelOutput = numRowSubArrayPhysical;                   // # of levels of the multilevelSenseAmp output, should be in 2^N forms; e.g. 32 levels --> 5-bit ADC
	cellBit = 1;                        // precision of memory device 
	// 1.4 update: dummy column sharing - how many senseamplfiers share one dummny columns?
	// dummy column sharing should not be high, since it could change the column cap of the dummy column.
	// In order for the dummy column to serve as a reference, the column caps of the dummy & main column should be matched closely
	dumcolshared = levelOutput;

	/*** parameters for SRAM ***/
	// due the scaling, suggested SRAM cell size above 22nm: 160F^2
	// SRAM cell size at 14nm: 300F^2
	// SRAM cell size at 10nm: 400F^2
	// SRAM cell size at 7nm: 600F^2
	heightInFeatureSizeSRAM = 10;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 28;        // SRAM Cell width in feature size  
	widthSRAMCellNMOS = 2;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;

	// 1.4 update : SRAM size update

	if (technode>14){
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 10;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 28;        // SRAM Cell width in feature size  
	}
	else if (technode==14){ // Samsung 14 nm 
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 10.6;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 30.8;        // SRAM Cell width in feature size  
	}
	else if (technode==10){ // TSMC 10 nm 
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 12.8;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 31.25;        // SRAM Cell width in feature size  
	}
	else if (technode==7){ // TSMC IEDM 2016
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 16;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 34.43;        // SRAM Cell width in feature size  
	}
	else if (technode==5){ // IRDS
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 19.2;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 43.75;        // SRAM Cell width in feature size  
	}
	else if (technode==3){ // IRDS
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 30;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 68.26;        // SRAM Cell width in feature size  
	}
	else if (technode==2){ // IRDS
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 42;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 120;// 111.42;        // SRAM Cell width in feature size  
	}
	else if (technode==1){ // IRDS
	widthSRAMCellNMOS = 1;                            
	widthSRAMCellPMOS = 1;
	widthAccessCMOS = 1;
	heightInFeatureSizeSRAM = 80;        // SRAM Cell height in feature size  
	widthInFeatureSizeSRAM = 144;        // SRAM Cell width in feature size  
	}


	minSenseVoltage = 0.1;
	
	/*** parameters for analog synaptic devices ***/
	heightInFeatureSize1T1R = 4;        // 1T1R Cell height in feature size
	widthInFeatureSize1T1R = 12;         // 1T1R Cell width in feature size
	heightInFeatureSizeCrossbar = 2;    // Crossbar Cell height in feature size
	widthInFeatureSizeCrossbar = 2;     // Crossbar Cell width in feature size
	
	resistanceOn = 6e3; // 6e3;               // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
	resistanceOff = 6e3*17;// 6e3*17;           // Roff resistance at Vr in the reported measurement dat (need to recalculate below if considering the nonlinearity)
	maxConductance = (double) 1/resistanceOn;
	minConductance = (double) 1/resistanceOff;
	
	// 230920 update 
	// read voltage needed for mux energy calculation - read voltage is fixed for Neurosim1.4, due to the assumptions for the ADC modeling equation (refer to manual for more information)
	if (technode == 130) {readVoltage=0.58;}
	else if (technode == 90) {readVoltage=0.58;}
	else if (technode == 65) {readVoltage=0.55;}
	else if (technode == 45) {readVoltage=0.51;}
	else if (technode == 32) {readVoltage=0.51;}
	else if (technode == 22) {readVoltage=0.51;}
	else if (technode == 14) {readVoltage=0.277;}
	else if (technode == 10) {readVoltage=0.28;} 
	else if (technode == 7) {readVoltage=0.264;}
	else if (technode == 5) {readVoltage=0.253;}
	else if (technode == 3) {readVoltage=0.248;}
	else if (technode == 2) {readVoltage=0.28;}
	else if (technode == 1) {readVoltage=0.272;}

	readPulseWidth = 10e-9;             // read pulse width in sec
	accessVoltage = 1.1;                // Gate voltage for the transistor in 1T1R
	resistanceAccess = resistanceOn*IR_DROP_TOLERANCE;            // resistance of access CMOS in 1T1R
	writeVoltage = 2;					// Enable level shifer if writeVoltage > 1.5V
	
	// 20250206 update
	if (memcelltype == 4) {		// for cap array
		heightInFeatureSize1T1R = 2;       // 1T1R Cell height in feature size
		widthInFeatureSize1T1R = 2;        // 1T1R Cell width in feature size
		resistanceOn = 1e20;               // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
		resistanceOff = 1e20*25;        // Roff resistance = Ron * onoff ratio
		maxConductance = (double) 1/resistanceOn;
		minConductance = (double) 1/resistanceOff;
		readVoltage = 0.06;	                // On-chip read voltage for memory cell
		accessVoltage = 0.06;                // Gate voltage for the transistor in 1T1R
		chargeDelay = 5e-9;				// Time to transfer charges to Cref
	}

	if (memcelltype == 5) {         // for 2TnC  array
                heightInFeatureSize2TnC = 4;       // 2TnC Cell height in feature size
                widthInFeatureSize2TnC = 12;        // 2TnC Cell width in feature size
                resistanceOn = 1e20;               // Ron resistance at Vr in the reported measurement data (need to recalculate below if considering the nonlinearity)
                resistanceOff = 1e20*25;        // Roff resistance = Ron * onoff ratio
                maxConductance = (double) 1/resistanceOn;
                minConductance = (double) 1/resistanceOff;
                accessVoltage = 1.2;                // Gate voltage for the transistor
                chargeDelay = 15e-9;     		// Time to transfer charges to Cref
		writePulseWidth = 10e-9;

		readDisturbFactor = 1e-6;   // 0.0001% shift per read
        	writeDisturbFactor = 1e-4;  // 0.01% shift per V/3 half-select

		{
			double vCoercive = ecFE * 1e8 * tFE;
			double vBiasNode = vthReadTr + readOverdrive;
			if      (mem_rdo == 2) readVoltage = 0.75*vCoercive + vBiasNode;  
			//else if (mem_rdo == 1) readVoltage = 0.25*vCoercive + vBiasNode;  
			else if (mem_rdo == 1) readVoltage = 0.75*vCoercive + vBiasNode;  
			else                   readVoltage = 0.75*vCoercive + vBiasNode;                
			// else                   readVoltage = writeVoltage;                
		}
		// if (mem_rdo == 1) {
		// 	readVoltage = 0.1;
		// } else if (mem_rdo == 2) {
		// 	readVoltage = 1.2;
		// } else if (mem_rdo == 3) {
		// 	readVoltage = 2;
		// }
        }

	if (memcelltype == 6) {         // for 1TnC array
    		    heightInFeatureSize1TnC = 4;       // 1TnC Cell height in feature size (pitch limited)
    		    widthInFeatureSize1TnC = 4;        // 1TnC Cell width in feature size
    		    resistanceOn = 1e20;               // Leakage resistance
    		    resistanceOff = 1e20*25;
    		    maxConductance = (double) 1/resistanceOn;
    		    minConductance = (double) 1/resistanceOff;
    		    readVoltage = 1.0;                 // Read voltage needs to be high enough to switch polarization
    		    accessVoltage = 1.2;               // Gate voltage for the single transistor
    		    chargeDelay = 15e-9;     		   // Time to dump charge onto BL
    		    writePulseWidth = 10e-9;


		{
			double vCoercive = ecFE * 1e8 * tFE;
			double vBiasNode = vthReadTr + readOverdrive;
			if      (mem_rdo == 2) readVoltage = 0.75*vCoercive + vBiasNode;  
			else if (mem_rdo == 1) readVoltage = 0.25*vCoercive + vBiasNode;  
			else                   readVoltage = writeVoltage;                
		}
		//     if (mem_rdo == 1) {
                // 	        readVoltage = 0.1;
                // 	} else if (mem_rdo == 2) {
                // 	        readVoltage = 1.2;
                // 	} else if (mem_rdo == 3) {
                // 	        readVoltage = 2;
                // 	}

    	}
	if (memcelltype == 7) {         // for 1T1C array
                    heightInFeatureSize1T1C = 4;       // 1T1C Cell height in feature size (pitch limited)
                    widthInFeatureSize1T1C = 4;        // 1T1C Cell width in feature size
                    resistanceOn = 1e20;               // Leakage resistance
                    resistanceOff = 1e20*25;
                    maxConductance = (double) 1/resistanceOn;
                    minConductance = (double) 1/resistanceOff;
                    readVoltage = 1.0;                 // Read voltage needs to be high enough to switch polarization
                    accessVoltage = 1.2;               // Gate voltage for the single transistor
                    chargeDelay = 15e-9;                    // Time to dump charge onto BL
                    writePulseWidth = 10e-9;

                    if (mem_rdo == 1) {
                                readVoltage = 0.1;
                        } else if (mem_rdo == 2) {
                                readVoltage = 1.2;
                        } else if (mem_rdo == 3) {
                                readVoltage = 2;
                        }

        }



	/*** Calibration parameters ***/
	if(validated){
		alpha = 1.44;	// wiring area of level shifter
		beta = 1.4;  	// latency factor of sensing cycle
		gamma = 0.5; 	// switching activity of DFF in shifter-add and accumulator
		delta = 0.15; 	// switching activity of adder 
		epsilon = 0.05; // switching activity of control circuits
		zeta = 1.22; 	// post-layout energy increase
	}		
	
	
	if (memcelltype == 1) {
		cellBit = 1;             // force cellBit = 8 for all SRAM cases
	} 
	

	
	/*** Initialize interconnect wires ***/

	// 1.4 update: wirewidth
	if (wireWidth >= 175) {
		AR = 1.6; 
		Rho = 2.01*1e-8;
	} else if ((110 <= wireWidth) &&  (wireWidth < 175)) {
		AR = 1.6; 
		Rho = 2.20*1e-8;
	} else if ((105 <= wireWidth) &&  (wireWidth < 110)) {
		AR = 1.7; 
		Rho = 2.21*1e-8;
	} else if ((80 <= wireWidth) &&  (wireWidth < 105)){
		AR = 1.7; 
		Rho = 2.37*1e-8;
	} else if ((56 <= wireWidth) &&  (wireWidth < 80)){
		AR = 1.8; 
		Rho = 2.63*1e-8;
	} else if ((40 <= wireWidth) &&  (wireWidth < 56)) {
		AR = 1.9; 
		Rho = 2.97*1e-8;
	} else if ((32 <= wireWidth) &&  (wireWidth < 40)) {
		AR = 2.0; 
		Rho = 3.25*1e-8;
	} else if ((22 <= wireWidth) &&  (wireWidth < 32)){
		AR = 2.00; Rho = 3.95*1e-8;
	} else if ((20 <= wireWidth) &&  (wireWidth < 22)){
		AR = 2.00; Rho = 4.17*1e-8; 
	} else if ((15 <= wireWidth) &&  (wireWidth < 20)){
		AR = 2.00; Rho = 4.98*1e-8; 
	} else if ((12 <= wireWidth) &&  (wireWidth < 15)){
		AR = 2.00; Rho = 5.8*1e-8; 
	} else if ((10 <= wireWidth) &&  (wireWidth < 12)){
		// AR = 3.00; Rho = 6.65*1e-8; 
		AR = 2.00; Rho = 6.61*1e-8; 
	} else if ((8 <= wireWidth) &&  (wireWidth < 10)){
		AR = 3.00; Rho = 7.87*1e-8; 
	} else {
		exit(-1); puts("Wire width out of range"); 
	}
	
	Rho = Rho * 1 / (1- ( (2*AR*wireWidth + wireWidth)*barrierthickness / (AR*pow(wireWidth,2) ) ));
	
	// 1.4 update: Metal0
	if (Metal0 >= 175) {
		AR_Metal0 = 1.6; 
		Rho_Metal0 = 2.01*1e-8;
	} else if ((110 <= Metal0) &&  (Metal0< 175)) {
		AR_Metal0 = 1.6; 
		Rho_Metal0 = 2.20*1e-8;
	} else if ((105 <= Metal0) &&  (Metal0< 110)){
		AR_Metal0 = 1.7; 
		Rho_Metal0 = 2.21*1e-8;
	} else if ((80 <= Metal0) &&  (Metal0< 105)) {
		AR_Metal0 = 1.7; 
		Rho_Metal0 = 2.37*1e-8;
	} else if ((56 <= Metal0) &&  (Metal0< 80)){
		AR_Metal0 = 1.8; 
		Rho_Metal0 = 2.63*1e-8;
	} else if ((40 <= Metal0) &&  (Metal0< 56)) {
		AR_Metal0 = 1.9; 
		Rho_Metal0 = 2.97*1e-8;
	} else if ((32 <= Metal0) &&  (Metal0< 40)) {
		AR_Metal0 = 2.0; 
		Rho_Metal0 = 3.25*1e-8;
	} else if ((22 <= Metal0) &&  (Metal0< 32)){
		AR_Metal0 = 2.00; Rho_Metal0 = 3.95*1e-8;
	} else if ((20 <= Metal0) &&  (Metal0< 22)){
		AR_Metal0 = 2.00; Rho_Metal0 = 4.17*1e-8; 
	} else if ((15 <= Metal0) &&  (Metal0< 20)){
		AR_Metal0 = 2.00; Rho_Metal0 = 4.98*1e-8; 
	} else if ((12 <= Metal0) &&  (Metal0< 15)){
		AR_Metal0 = 2.00; Rho_Metal0 = 5.8*1e-8; 
	} else if ((10 <= Metal0) &&  (Metal0< 12)){
		// AR_Metal0 = 3.00; Rho_Metal0 = 6.65*1e-8; 
		AR_Metal0 = 2.00; Rho_Metal0 = 6.61*1e-8; 
	} else if ((8 <= Metal0) &&  (Metal0< 10)){
		AR_Metal0 = 3.00; Rho_Metal0 = 7.87*1e-8; 
	} else {
		exit(-1); puts("Wire width out of range"); 
	}
	
	Rho_Metal0 = Rho_Metal0 * 1 / (1- ( (2*AR_Metal0*Metal0 + Metal0)*barrierthickness / (AR_Metal0*pow(Metal0,2) ) ));
	
	
	// 1.4 update: Metal1
	if (Metal1 >= 175) {
		AR_Metal1 = 1.6; 
		Rho_Metal1 = 2.01*1e-8;
	} else if ((110 <= Metal1) &&  (Metal1 < 175)) {
		AR_Metal1 = 1.6; 
		Rho_Metal1 = 2.20*1e-8;
	} else if ((105 <= Metal1) &&  (Metal1 < 110)) {
		AR_Metal1 = 1.7; 
		Rho_Metal1 = 2.21*1e-8;
	} else if ((80 <= Metal1) &&  (Metal1 <105)) {
		AR_Metal1 = 1.7; 
		Rho_Metal1 = 2.37*1e-8;
	} else if ((56 <= Metal1) &&  (Metal1 < 80)) {
		AR_Metal1 = 1.8; 
		Rho_Metal1 = 2.63*1e-8;
	} else if ((40 <= Metal1) &&  (Metal1 < 56)){
		AR_Metal1 = 1.9; 
		Rho_Metal1 = 2.97*1e-8;
	} else if ((32 <= Metal1) &&  (Metal1 < 40)) {
		AR_Metal1 = 2.0; 
		Rho_Metal1 = 3.25*1e-8;
	} else if ((22 <= Metal1) &&  (Metal1 < 32)){
		AR_Metal1 = 2.00; Rho_Metal1 = 3.95*1e-8;
	} else if ((20 <= Metal1) &&  (Metal1 < 22)){
		AR_Metal1 = 2.00; Rho_Metal1 = 4.17*1e-8; 
	} else if ((15 <= Metal1) &&  (Metal1 < 20)){
		AR_Metal1 = 2.00; Rho_Metal1 = 4.98*1e-8; 
	} else if ((12 <= Metal1) &&  (Metal1 < 15)){
		AR_Metal1 = 2.00; Rho_Metal1 = 5.8*1e-8; 
	} else if ((10 <= Metal1) &&  (Metal1 < 12)){
		// AR_Metal1 = 3.00; Rho_Metal1 = 6.65*1e-8; 
		AR_Metal1 = 2.00; Rho_Metal1 = 6.61*1e-8;
	} else if ((8 <= Metal1) &&  (Metal1 < 10)){
		AR_Metal1 = 3.00; Rho_Metal1 = 7.87*1e-8; 
	} else {
		exit(-1); puts("Wire width out of range"); 
	}

	Rho_Metal1 = Rho_Metal1 * 1 / (1- ( (2*AR_Metal1*Metal1 + Metal1)*barrierthickness / (AR_Metal1*pow(Metal1,2) ) ));
	
	Metal0_unitwireresis =  Rho_Metal0 / ( Metal0*1e-9 * Metal0*1e-9 * AR_Metal0 );
	Metal1_unitwireresis =  Rho_Metal1 / ( Metal1*1e-9 * Metal1*1e-9 * AR_Metal1 );
	


	if (memcelltype == 1) {
		wireLengthRow = wireWidth * 1e-9 * heightInFeatureSizeSRAM;
		wireLengthCol = wireWidth * 1e-9 * widthInFeatureSizeSRAM;
	} else {
		if (accesstype == 1) {
			wireLengthRow = wireWidth * 1e-9 * heightInFeatureSize1T1R;
			wireLengthCol = wireWidth * 1e-9 * widthInFeatureSize1T1R;
		} else {
			wireLengthRow = wireWidth * 1e-9 * heightInFeatureSizeCrossbar;
			wireLengthCol = wireWidth * 1e-9 * widthInFeatureSizeCrossbar;
		}
	}
	Rho *= (1+0.00451*abs(temp-300));
	if (wireWidth == -1) {
		unitLengthWireResistance = 1.0;	// Use a small number to prevent numerical error for NeuroSim
		wireResistanceRow = 0;
		wireResistanceCol = 0;
	} else {
		unitLengthWireResistance =  Rho / ( wireWidth*1e-9 * wireWidth*1e-9 * AR );
		wireResistanceRow = unitLengthWireResistance * wireLengthRow;
		wireResistanceCol = unitLengthWireResistance * wireLengthCol;
	}
}

