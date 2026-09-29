/******************out*************************************************************
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

#ifndef PARAM_H_
#define PARAM_H_

class Param {
public:
	Param();

	int operationmode, operationmodeBack, memcelltype, accesstype, transistortype, deviceroadmap;     
	
	int mem_rdo;	
	
	double heightInFeatureSizeSRAM, widthInFeatureSizeSRAM, widthSRAMCellNMOS, widthSRAMCellPMOS, widthAccessCMOS, minSenseVoltage;
	 
	double heightInFeatureSize1T1R, widthInFeatureSize1T1R, heightInFeatureSizeCrossbar, widthInFeatureSizeCrossbar;
	
	int relaxArrayCellHeight, relaxArrayCellWidth;
	// Anni update
	bool globalBusType, globalBufferType, tileBufferType, peBufferType, chipActivation, reLu, novelMapping, pipeline, SARADC, currentMode, validated, synchronous;
	int globalBufferCoreSizeRow, globalBufferCoreSizeCol, tileBufferCoreSizeRow, tileBufferCoreSizeCol;																								
	
	double clkFreq, featuresize, readNoise, resistanceOn, resistanceOff, maxConductance, minConductance;
	int temp, technode, wireWidth, multipleCells;
	double maxNumLevelLTP, maxNumLevelLTD, readVoltage, readPulseWidth, writeVoltage;
	double accessVoltage, resistanceAccess;
	double nonlinearIV, nonlinearity;
	double writePulseWidth, numWritePulse;
	double globalBusDelayTolerance, localBusDelayTolerance;
	double treeFoldedRatio, maxGlobalBusWidth;
	double algoWeightMax, algoWeightMin;
	
	int neuro, multifunctional, parallelWrite, parallelRead;
	int numlut, numColMuxed, numWriteColMuxed, levelOutput, avgWeightBit, numBitInput;
	int numRowSubArray, numColSubArray;
	int cellBit, synapseBit;
	int speedUpDegree;
	
	int XNORparallelMode, XNORsequentialMode, BNNparallelMode, BNNsequentialMode, conventionalParallel, conventionalSequential; 
	int numRowPerSynapse, numColPerSynapse;
	double AR, Rho, wireLengthRow, wireLengthCol, unitLengthWireResistance, wireResistanceRow, wireResistanceCol;
	
	double alpha, beta, gamma, delta, epsilon, zeta;

	// 1.4 update: BEOL related parameters added
	
	double Metal0=0;
	double Metal1=0;
	double AR_Metal0=0;
	double AR_Metal1=0;
	double Rho_Metal0=0;
	double Rho_Metal1=0;
	double Metal0_unitwireresis=0;
	double Metal1_unitwireresis=0;

	// 1.4 update: add activation implementation option

	bool Activationtype; // true: SRAM, False: RRAM

	// 1.4 update: Final driver sizing for row decoder conventional parallel mode (SRAM, RRAM)
	// multiplied by the driver width
	double sizingfactor_MUX= 1; 
	double sizingfactor_WLdecoder= 1; 

	// 1.4 update: switchmatrix parameter tuning
	double newswitchmatrixsizeratio=6;
	double switchmatrixsizeratio=1;
	
	// 1.4 update: Special layout
	double speciallayout;
	
	// 1.4 update: added parameters for buffer insertion
	double unitcap;
	double unitres;
	double drivecapin; 
	double buffernumber=0;
	double buffersizeratio=0;
	
	// 1.4 update: barrier thickness
	double barrierthickness= 0;
	
	// 1.4 update: new ADC modeling related parameters
	double dumcolshared;
	double columncap;
	double reference_energy_peri=0;
	
	// 1.4 update: array dimension/SRAM access resistance for multilevelsenseamp
	double arrayheight;
	double arraywidthunit;
	double resCellAccess;

	// 1.4 update 
	double inputtoggle;
	double outputtoggle;

	// 1.4 debug
	double ADClatency;
	double rowdelay;
	double muxdelay;
	
	// 1.4 update: technology node
	int technologynode;

	// Anni update: partial parallel mode
	int numRowParallel;

	// 230920 update
	double totaltile_num;
	int sync_data_transfer;

	// Cap update 20250206
	double chargeDelay;


	double heightInFeatureSize2TnC;
        double widthInFeatureSize2TnC;
        
	double heightInFeatureSize1TnC;
        double widthInFeatureSize1TnC;

	double heightInFeatureSize1T1C;
        double widthInFeatureSize1T1C;

	int numCapacitors;
        double capacitance;

	int numRowSubArrayPhysical;
	int numColSubArrayPhysical;
	int bitsPerCell; // Set this to 8 for your 2T8C
	
	double readDisturbFactor;
    	double writeDisturbFactor;

	// 3D floorplan: architecture + staircase 
	int    integrationMode;        // 0=CNA, 1=CUA, 2=CBA
	double staircaseStepWidth;     // minimum stair-tread run                          [m]
	double viaDiameter3D;          // plane/WL contact via diameter                    [m]
	double viaOverlayMargin;       // landing-pad overlay margin per side              [m]
	double viaSpacing3D;           // edge-to-edge keep-out between via pads           [m]
	double metalPitch3D;           // BEOL metal pitch for per-layer routing           [m]
	double logicPackingEfficiency; // logic-die fill factor (0-1)
	double bondPadPitch;           // Cu-Cu hybrid-bond pad pitch (CBA)                [m]
	int    numBondPadPerSubarray;  // die-to-die bond pads per subarray (CBA)

	double staircaseStepPitch;
	int    staircaseDummySteps;
	bool   staircaseBothSidesFullContact;
	double staircaseEdgeMargin;
	double tavPitch;
	double cuaUtilization;
	double cuaAreaDerate;
	double cbaUtilization;

	// memory operator mode 
	bool   memoryMode;              // true: pure memory macro, no CIM dataflow
	bool   coreOnly;                // true: also accumulate core-array-only PPA
	int    numRowActivated;         // word lines asserted per access (1 = memory)
	double qndroRefreshInterval;    // N: qndro reads tolerated before a page rewrite

	/* ---- sense tier: core array + column sense amplifier, no row drivers ---- */
	double senseLatencyCore;      // one sense event, worst column
	double senseEnergyCore;       // all numCol columns, one access
	double readLatencySense;      // = readLatencyCore + senseLatencyCore
	double readEnergySense;       // = readEnergyCore  + senseEnergyCore

	double twoPr;          // remanent polarisation 2Pr, C/m^2
	double cellAreaFE;     // ferroelectric capacitor electrode area, m^2

		/* ---- ferroelectric material, from measurement, not guesses ---- */
	double epsFE;              // relative permittivity of the HZO
	double tFE;                // ferroelectric thickness (m)
	/* ---- read transistor ---- */
	double vthReadTr;          // threshold voltage of TR (V)
	double ssReadTr;           // subthreshold swing (V/decade)
	double ionSatReadTr;       // saturated ON current (A)
	double ioffReadTr;         // OFF current (A)
	double capGateRatio;       // C_MFM / C_gate design ratio
	double senseMarginVolt;    // required |Vint - Vth| margin (V)

	/* ---- NLS polarisation switching kinetics ----
	 * tau(E) = nlsTauInf * exp[(nlsEa/E)^nlsAlpha],  E in MV/cm
	 * Alessandri et al., IEEE EDL 2018, 8 nm Hf0.5Zr0.5O2            [4] */
	double nlsTauInf;          // s      asymptotic switching time (HARD floor)
	double nlsEa;              // MV/cm  activation field
	double nlsAlpha;           // -      exponent
	double ecFE;               // MV/cm  coercive field, for the disturb check

	/* ---- read bias and sensing ---- */
	double readOverdrive;      // Vint_bias = vthReadTr + this (V)
	double senseAmpResolution; // differential input the sense amp resolves (V)

	/* ---- data pattern ----
	 * fraction of the accessed page holding the SWITCHING state.
	 * 0.5 reproduces the old hardcoded 50/50; 0.0 and 1.0 are the bounds. */
	double dataOnesRead;

	double inhibitDivider;       // 2.0 = V/2 self-boosted, 3.0 = V/3
	double qndroSwitchFraction;  // fraction of 2Pr switched by a QNDRO read
	double plateSheetRes;   // ohm/square of the WBL plate metal

	double resPlateDriver;
	double staircaseStep;
	double staircasePadPitch;

	int    inhibitPillarModel;   /* 0 = unselected pillars biased at Vw/2 (V/2 scheme)
	                              * 1 = unselected pillars float (self-boost, pessimistic) */
	bool   inhibitSupplyDraw;    /* true = C*V^2 (supply), false = 0.5*C*V^2 (stored) */

		/* ---- read-transistor bias and sense amplifier ---- */
	bool   readOverdriveAuto;  /* true: bias = window - dV/2, set per run in rw_main */
	double senseAmpVmin;       /* differential the LATCH needs at its own input (V) */
	double capSenseInt;        /* dedicated integrating cap; 0 => integrate on capRBL */

	int    senseMode;          /* 0 = clamped column, current-mode; 1 = column integrates */
	double iClampBias;         /* clamp bias current (A): g_clamp = iClampBias/(n*kT/q) */
	double capJunctionTr;      /* S/D capacitance one transistor adds to its line (F); 0 = NeuroSim planar model */
	double capSenseAmpIn;      /* amplifier input capacitance added to the integrating node (F) */
	double capSenseAmpOut;     /* load on the latch's regenerating node: output register + local wiring (F) */
	double senseRefFraction;   /* fraction of (I1-I0) available as differential: 0.5 = mid-point reference */
	double senseAmpRegenTime;  /* latch regeneration time (s); 0 = derive tau = C/gm from the technology */
	bool   readDrivesUnselPlanes; /* true: the (n-1) unselected strips of the accessed row are driven to vUnsel per read */
	double kLineRc;            /* driven-line settle criterion in units of R*C (1.0 distributed 90 %, 2.2 lumped) */
	double rblResScale;        /* RBL resistance relative to the technology minimum wire (0.25 = 4x wider or double-ended) */

};

#endif
