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
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>
#include "constant.h"
#include "formula.h"
#include "SubArray.h"
#include "BlockStats.h"
#include "Param.h"
//#include "Floorplan3D.h"


using namespace std;


static const double EPS0_VACUUM = 8.854e-12;   // F/m

extern Param *param;

SubArray::SubArray(InputParameter& _inputParameter, Technology& _tech, MemCell& _cell):
						inputParameter(_inputParameter), tech(_tech), cell(_cell),
						wllevelshifter(_inputParameter, _tech, _cell),
						sllevelshifter(_inputParameter, _tech, _cell),
						bllevelshifter(_inputParameter, _tech, _cell),
						blDecoder(_inputParameter, _tech, _cell),
						plDecoder(_inputParameter, _tech, _cell),
						wlDecoder(_inputParameter, _tech, _cell),
						wlDecoderOutput(_inputParameter, _tech, _cell),
						wlNewDecoderDriver(_inputParameter, _tech, _cell),
						wlNewSwitchMatrix(_inputParameter, _tech, _cell),
						rowCurrentSenseAmp(_inputParameter, _tech, _cell),
						mux(_inputParameter, _tech, _cell),
						muxDecoder(_inputParameter, _tech, _cell),
						slSwitchMatrix(_inputParameter, _tech, _cell),
						blSwitchMatrix(_inputParameter, _tech, _cell),
						wlSwitchMatrix(_inputParameter, _tech, _cell),
						deMux(_inputParameter, _tech, _cell),
						readCircuit(_inputParameter, _tech, _cell),
						precharger(_inputParameter, _tech, _cell),
						senseAmp(_inputParameter, _tech, _cell),
						wlDecoderDriver(_inputParameter, _tech, _cell),
						sramWriteDriver(_inputParameter, _tech, _cell),
						adder(_inputParameter, _tech, _cell),
						dff(_inputParameter, _tech, _cell),
						shiftAddInput(_inputParameter, _tech, _cell),
						shiftAddWeight(_inputParameter, _tech, _cell),
						currentSenseAmp(_inputParameter, _tech, _cell),
						multilevelSenseAmp(_inputParameter, _tech, _cell),
						multilevelSAEncoder(_inputParameter, _tech, _cell),
						plSwitchMatrix(_inputParameter, _tech, _cell),
						wplSwitchMatrix(_inputParameter, _tech, _cell),
						rslSwitchMatrix(_inputParameter, _tech, _cell),
						rblSwitchMatrix(_inputParameter, _tech, _cell),
						sslSwitchMatrix(_inputParameter, _tech, _cell),
                                                wblSwitchMatrix(_inputParameter, _tech, _cell),
                                                wwlSwitchMatrix(_inputParameter, _tech, _cell),
						wwlDecoder(_inputParameter, _tech, _cell),
						wplDecoder(_inputParameter, _tech, _cell),
						wblDecoder(_inputParameter, _tech, _cell),
						wblPlaneDecoder(_inputParameter, _tech, _cell),
						rslDecoder(_inputParameter, _tech, _cell),
						rblDecoder(_inputParameter, _tech, _cell),
						sslDecoder(_inputParameter, _tech, _cell),
						blLevelShifter(_inputParameter, _tech, _cell),
						wblLevelShifter(_inputParameter, _tech, _cell),
                        			sslLevelShifter(_inputParameter, _tech, _cell),
                        			wlLevelShifter(_inputParameter, _tech, _cell),
                        			wwlLevelShifter(_inputParameter, _tech, _cell),
                        			plLevelShifter(_inputParameter, _tech, _cell),
                        			wplLevelShifter(_inputParameter, _tech, _cell),
						sarADC(_inputParameter, _tech, _cell){
						
						initialized = false;
						readDynamicEnergyArray = writeDynamicEnergyArray = 0;
						activityRowRead  = 0;          
						activityRowWrite = 0;          
						activityColRead  = 1;          
						activityColWrite = 1;
						readLatencyCore = writeLatencyCore = 0;
						lineSetupCore = restoreLatencyCore = 0;
						readEnergyCore = readSenseEnergyCore = readRestoreEnergyCore = 0;
						writeEnergyCore = writeCellEnergyCore = inhibitionEnergyCore = 0;
						areaCore = 0;
						senseLatencyCore = senseEnergyCore = 0;
						readLatencySense = readEnergySense = 0;

						gateCap = 0; 
						capWBLwire = capWBLpar = 0; 
						capRSL = capWBL = capRBL = capSSL = capWPL = capWWL = 0;	
} 

void SubArray::Initialize(int _numRow, int _numCol, double _unitWireRes){  //initialization module
	
	//numRow = _numRow;    // import parameters
	
	if (cell.memCellType == Type::_2TnC || cell.memCellType == Type::_1TnC) {
		numRow = param->numRowSubArrayPhysical; 
	} else {
		numRow = _numRow;    // import parameters
	}

	numCol = _numCol;
	unitWireRes = _unitWireRes;
	double bitsPerCell = param->bitsPerCell;
	
	double MIN_CELL_HEIGHT = MAX_TRANSISTOR_HEIGHT;  //set real layout cell height
	double MIN_CELL_WIDTH = (MIN_GAP_BET_GATE_POLY + POLY_WIDTH) * 2;  //set real layout cell width
	double ISOLATION_REGION = MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY; // 1.4 update : new variable

	// 1.4 update : new cell dimension setting

	if (tech.featureSize == 14 * 1e-9){
	MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_14nm/MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_14nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 10 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_10nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_10nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 7 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_7nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_7nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 5 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_5nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_5nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 3 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_3nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_3nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 2 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_2nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_2nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else if (tech.featureSize == 1 * 1e-9){
    MIN_CELL_HEIGHT *= (MAX_TRANSISTOR_HEIGHT_1nm /MAX_TRANSISTOR_HEIGHT);
	ISOLATION_REGION *= (OUTER_HEIGHT_REGION_1nm/(MIN_POLY_EXT_DIFF*2 + MIN_GAP_BET_FIELD_POLY));}
    else{
    MIN_CELL_HEIGHT *= 1;
	ISOLATION_REGION *=1;}

	if (tech.featureSize == 14 * 1e-9)
	MIN_CELL_WIDTH  *= ((POLY_WIDTH_FINFET + MIN_GAP_BET_GATE_POLY_FINFET )/(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 10 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_10nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 7 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_7nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 5 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_5nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 3 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_3nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 2 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_2nm /(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else if (tech.featureSize == 1 * 1e-9)
    MIN_CELL_WIDTH  *= (CPP_1nm/(MIN_GAP_BET_GATE_POLY + POLY_WIDTH));
    else
    MIN_CELL_WIDTH  *= 1;

	// 1.4 update : think about the factor multiplied by the cell width/height for the relaxation 

	if (cell.memCellType == Type::SRAM) {  //if array is SRAM
		if (relaxArrayCellWidth) {  //if want to relax the cell width

			// 1.4 update: For SRAM, the height corresponds to the width of the logic layout 

			lengthRow = (double)numCol * MAX(cell.widthInFeatureSize, MIN_CELL_HEIGHT) * tech.featureSize;
		} else { //if not relax the cell width
			lengthRow = (double)numCol * cell.widthInFeatureSize * tech.featureSize;
		}
		if (relaxArrayCellHeight) {  //if want to relax the cell height

			// 1.4 update: For SRAM, the height corresponds to the width of the logic layout 

			lengthCol = (double)numRow * MAX(cell.heightInFeatureSize, MIN_CELL_WIDTH) * tech.featureSize;
		} else {  //if not relax the cell height
			lengthCol = (double)numRow * cell.heightInFeatureSize * tech.featureSize;
		}

		// 230920 update
		param->arraywidthunit = cell.widthInFeatureSize * tech.featureSize;
		param->arrayheight = (double)numRow * cell.heightInFeatureSize * cell.featureSize;


	} else if (cell.memCellType == Type::_2TnC) {
                // 2TnC Initialization
                // double cellHeight = cell.heightInFeatureSize;
                // double cellWidth = cell.widthInFeatureSize;

                //Calculate Array Dimensions
                // if (relaxArrayCellWidth) {
                //      lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;
                // } else {
                //      lengthRow = (double)numCol * cellWidth * tech.featureSize;
                // }

                // if (relaxArrayCellHeight) {
                //      lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
                // } else {
                //      lengthCol = (double)numRow * cellHeight * tech.featureSize;
                // }


                // Calculate Parasitics
                // double gateCap = CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech);
		
		// Cgate = (2 * pi * epsilon_ox * L_g) / ln(1 + (tox/r)) 
		// r_wire = 4e-9 // 4nm nanowire radius
		// t_ox = 1.2e-9 // 1.2nm oxide thickness
		// L_g = 22e-9   // 22nm gate length
		// eps_ox = 3.9 * 8.85e-12 // SiO2 permittivity 

                //double gateCap = 0.017e-15;
                double gateCap = 0.026e-15;
                // double drainCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);
		// double sourceCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);

                this->gateCap = gateCap;
                double drainCap  = (param->capJunctionTr > 0.0) ? param->capJunctionTr
                                 : CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);
		double sourceCap = (param->capJunctionTr > 0.0) ? param->capJunctionTr
		                 : CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);

                // 3D 2T-nC WBL GEOMETRY CALCULATION

                // 1. Define Dimensions (in m)
                double tFE = 5e-9;          // FE thickness
                double rString = 40e-9;     // String radius (XY)
                double distString = 40e-9;  // String interval (XY)
                double tWBL = 40e-9;        // WBL thickness (Z)
                double distWBL = 20e-9;     // WBL interval (Z) - Spacing between layers

                // Calculate Pitch (Unit Cell Size)
                // Pitch = Diameter of String + Spacing (Footprint per cell)
                double pitch3D = (2 * rString) + distString;

		// 2. Override Array Dimensions for 3D
        	// In 3D Vertical arrays, the "Cell Width/Height" is the pillar pitch
        	// We ignore relaxArrayCellWidth for 3D as it's pitch-limited
        	lengthRow = (double)numCol * pitch3D;
        	lengthCol = (double)numRow * pitch3D;

		// Update global params for the array size
        	param->arraywidthunit = pitch3D;
        	param->arrayheight = lengthCol;

		// 3. WBL Calculations
		// 3D WBL Resistance
        	double rho = 10e-8;
        	double holeDiameter = 2 * (rString + tFE);  //Hole Diameter includes the FE layer thickness around the string 
        	double effectiveWidth = pitch3D - holeDiameter; //Effective Width = Pitch - Hole Diameter (The narrowest path for current)
        	double resPerCellWBL = rho * (pitch3D / (effectiveWidth * tWBL));  //Resistance = Rho * (Length / Area_cross_section)

        	// 3D WBL Capacitance  (Parasitic + Cell)
        	double areaUnitCell = pitch3D * pitch3D;
        	double areaHole = 3.14159 * pow(rString + tFE, 2);
        	double areaMetal = areaUnitCell - areaHole;  //Area of Metal Plane = (Pitch^2) - (Area of Hole)
        	double epsilonOxide = 3.9 * 8.85e-12;
        	double capCouplingPerCell = epsilonOxide * areaMetal / distWBL;
        	double capParasiticWBL = 2 * capCouplingPerCell * numCol * numRow; // Multiply by 2 because it couples to BOTH the layer above and below
										   // capParasiticWBL is the single cell capacitance, so we multiply it with numCol * numRow

		double epsilonFE = 25 * 8.85e-12;
        	double capCell3D = (2 * 3.14159 * epsilonFE * tWBL) / log((rString + tFE) / rString);
		double nStack   = (double)bitsPerCell;                               // planes per pillar
        	double capSer   = (nStack > 1.0)
        	                ? capCell3D * ((nStack - 1.0) * capCell3D + gateCap)
        	                            / ( nStack        * capCell3D + gateCap)
        	                : capCell3D;
										   
										   // Cell Capacitance (FeCAP)
        	double nSer = (bitsPerCell > 1)
	        	    ? ((double)bitsPerCell - 1.0) / (double)bitsPerCell : 1.0;
		///double capWBLTotal = capParasiticWBL + (capCell3D * nSer * numRow * numCol);
		//double capWBLTotal = capParasiticWBL + (capCell3D * numCol);

		double capWBLTotal = capParasiticWBL + (capSer * numCol * numRow);

                //double resRowWBL = resPerCellWBL * numCol;

		// 4. 3D Via & Crossing Caps 
        	// Defined trans_height as 2D feature size approx
        	double trans_height = 2 * MAX_TRANSISTOR_HEIGHT * tech.featureSize; // Write and Read Transistors
        	double distSSL = 200e-9;
		double distRSL = 100e-9;
		double viaHeight = ((tWBL + distWBL) * bitsPerCell) + distSSL + distRSL + trans_height;  // Height: Must span all 8 layers + spacing
        	double viaRadius = 40e-9;  // Same as String Radius

        	// RBL Crossing (Top)
        	double areaOverlap = (tech.featureSize * tech.featureSize);
        	double capCrossingRBL = epsilonOxide * areaOverlap / viaHeight;
        	double totalCapCrossRBL = (3 * capCrossingRBL);

        	// RSL Crossing (Top)
        	double capCrossingRSL = epsilonOxide * areaOverlap / viaHeight;
        	double totalCapCrossRSL = (3 * capCrossingRSL);

        	// WBL Coupling (Middle)
        	double h_coupling = tWBL;
        	double r_inner = rString;
        	double r_outer = rString + tFE;
        	double capCouplingPerPillar = (2 * 3.14159 * epsilonOxide * h_coupling) / log(r_outer / r_inner);
        	double totalCapCrossWBL = (2 * capCouplingPerPillar) * numRow * numCol;
		//double totalCapCrossWBL = (2 * capCouplingPerPillar) * numCol;


        	double lenStair  = (double)bitsPerCell * param->staircaseStepPitch;
        	double densStair = 2.0 * epsilonOxide / distWBL;      // F/m^2, both neighbours 
        	double capStair  = lenStair * lengthRow * densStair;

        	double sides     = 2.0;                               // staircase on 2 edges 
        	double sqArray   = lengthCol / lengthRow;             // squares across the array 
        	double sqStair   = lenStair  / lengthRow;             // squares along the tail   

        	capWBL = capWBLTotal + totalCapCrossWBL + capStair;
        	resWBL = param->plateSheetRes * (sqArray + sqStair) / sides;


        	// Summing total crossing cap for general line usage
        	double totalCrossingCapPerCol = totalCapCrossRBL + totalCapCrossRSL + totalCapCrossWBL;

        	// 5. Final Capacitance Assignments
		// WPL (Source Connected) - capCol
		capWPL = lengthRow * 0.2e-15/1e-6 + (sourceCap * numCol);

        	// WWL (Gate connected) - capColGate
        	capWWL = lengthRow * 0.2e-15/1e-6 + (gateCap * numCol);

        	// SSL (Source Select) - capRow1
        	capSSL = lengthCol * 0.2e-15/1e-6 + (sourceCap * numRow);

        	// WBL (Plates)
		double capPlateToPlate = 2.0 * epsilonOxide * (lengthRow * lengthCol) / distWBL;
		//capWBL   = capPlateToPlate + capWBLTotal + totalCapCrossWBL;
		//resWBL   = param->plateSheetRes;     
		//resWBL   = resRow;     
		//capWBL = lengthRow * 0.2e-15/1e-6 + (capWBLTotal + totalCapCrossWBL);

		capWBL = capWBLTotal + totalCapCrossWBL;

		double plateSquares = (double)numRow / (double)numCol;
        	resWBL = param->plateSheetRes * plateSquares
        	       / (param->plateContactSides > 0 ? param->plateContactSides : 1.0);

        	capWBLwire = lengthRow * 0.2e-15/1e-6;
		capWBLpar  = capWBLwire + capParasiticWBL + totalCapCrossWBL;   // capWBL minus numCol*capCell3D 

		// RSL (Top Select)
        	double wireCapCol = lengthCol * 0.2e-15/1e-6;
        	capRSL = wireCapCol + ((sourceCap + totalCapCrossRSL) * numRow);
        	
		// RBL 
		double DrainCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, pitch3D, tech);
		capRBL = lengthRow * 0.2e-15/1e-6 + ((DrainCap + totalCapCrossRBL) * numCol);

        	resRow = lengthRow * unitWireRes;
        	resCol = lengthCol * unitWireRes;


		double unitcap= capRBL/param->numColSubArray;
                double unitres= resRow/param->numColSubArray;
                param->unitcap = unitcap;
                param->unitres = unitres;

		param->columncap = capRSL;


                // Define Transmission Gate Resistance for Switch Matrices
                //double resTg = cell.resMemCellOn;

	
	
	} else if (cell.memCellType == Type::_1TnC) {
        	double gateCap = 0.017e-15;
                double drainCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);
                double sourceCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);

                // 3D 2T-nC WBL GEOMETRY CALCULATION
                // 1. Define Dimensions (in m)
                double tFE = 5e-9;          // FE thickness
                double rString = 40e-9;     // String radius (XY)
                double distString = 40e-9;  // String interval (XY)
                double tWBL = 40e-9;        // WBL thickness (Z)
                double distWBL = 20e-9;     // WBL interval (Z) - Spacing between layers

                // Calculate Pitch (Unit Cell Size)
                // Pitch = Diameter of String + Spacing (Footprint per cell)
                double pitch3D = (2 * rString) + distString;

                // 2. Override Array Dimensions for 3D
                // In 3D Vertical arrays, the "Cell Width/Height" is the pillar pitch
                // We ignore relaxArrayCellWidth for 3D as it's pitch-limited
                lengthRow = (double)numCol * pitch3D;
                lengthCol = (double)numRow * pitch3D;

                // Update global params for the array size
                param->arraywidthunit = pitch3D;
                param->arrayheight = lengthCol;

                // 3. WBL Calculations
                // 3D WBL Resistance
                double rho = 10e-8;
                double holeDiameter = 2 * (rString + tFE);  //Hole Diameter includes the FE layer thickness around the string
                double effectiveWidth = pitch3D - holeDiameter; //Effective Width = Pitch - Hole Diameter (The narrowest path for current)
                double resPerCellWBL = rho * (pitch3D / (effectiveWidth * tWBL));  //Resistance = Rho * (Length / Area_cross_section)

                // 3D WBL Capacitance  (Parasitic + Cell)
                double areaUnitCell = pitch3D * pitch3D;
                double areaHole = 3.14159 * pow(rString + tFE, 2);
                double areaMetal = areaUnitCell - areaHole;  //Area of Metal Plane = (Pitch^2) - (Area of Hole)
                double epsilonOxide = 3.9 * 8.85e-12;
                double capCouplingPerCell = epsilonOxide * areaMetal / distWBL;
                double capParasiticWBL = 2 * capCouplingPerCell * numCol; //Multiply by 2 because it couples to BOTH the layer above and below

                // Cell Capacitance (FeCAP)
                double epsilonFE = 25 * 8.85e-12;
                double capCell3D = (2 * 3.14159 * epsilonFE * tWBL) / log((rString + tFE) / rString);
                double capWBLTotal = capParasiticWBL + (capCell3D * numCol);
                // double capWBLTotal = capParasiticWBL;


                //double resRowWBL = resPerCellWBL * numCol;

		// PL 
                capPL = lengthCol * 0.2e-15/1e-6 + (sourceCap * numRow);

                // WL 
                capWL = lengthRow * 0.2e-15/1e-6 + (gateCap * numCol);

                // BL (Plates)
                capBL = lengthRow * 0.2e-15/1e-6 + (capWBLTotal);

        	resRow = lengthRow * unitWireRes;
        	resCol = lengthCol * unitWireRes;

		double unitcap= capPL/param->numColSubArray;
                double unitres= resRow/param->numColSubArray;
                param->unitcap = unitcap;
                param->unitres = unitres;

                param->columncap = capPL;


        	if (cell.mem_rdo == Type::ndro) {

                           // 1. Calculate the capacitive load at the internal node (Gate of Tr)
                           //double cGateTr = CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech);
                           double cGateTr = 0.017e-15;

                           // Define FeCap parameters
                           // These represent the switching capacitance differences.
                           double capFeOn = 10e-15;  // High switching capacitance (State '0')
                           double capFeOff = 2e-15;  // Low switching capacitance (State '1')

                           // 2. Voltage Coupling: Calculate V_int (Gate Voltage) for both states
                           double vGateOn = cell.readVoltage * (capFeOn / (capFeOn + cGateTr));
                           double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));

                           // 3. Extract Transistor parameters for physics equations

                           //Add custom values, tech is not the same
                           //double muCox = tech.mobility * tech.cox;
                           double muCox = 150e-6;
                           double w_l = cell.widthAccessCMOS;
                           //double vTh = tech.vth;
                           double vTh = 0.3;
                           double vT = 0.026; // Thermal voltage at room temp (kT/q)

                           // 4. Calculate "On" State Resistance (Linear Region of Tr)
                           if (vGateOn > vTh) {
                               cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
                           } else {
                               cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
                           }

                           // 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
                           double n_factor = 1.2; // Subthreshold swing factor
                           double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));

                           cell.resMemCellOff = cell.readVoltage / iSubOff;

		} else if (cell.mem_rdo == Type::qndro) {

                           // 1. In QNDR, the read voltage is high enough to flip the domain, dumping
                           // the full polarization charge onto the gate. This is represented with a massive capFeOn.

                           double cGateTr = 0.017e-15;

                           // QNDR Parameters
                           double capFeOn_QNDR = 40e-15;  // Massive effective switching capacitance
                           double capFeOff = 2e-15;       // Unswitched state (dielectric only)

                           // 2. Calculate Voltage Divider
                           double vGateOn = cell.readVoltage * (capFeOn_QNDR / (capFeOn_QNDR + cGateTr));
                           double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));

                           // 3. Extract Transistor parameters for physics equations

                           //Add custom values, tech is not the same
                           //double muCox = tech.mobility * tech.cox;
                           double muCox = 150e-6;
                           double w_l = cell.widthAccessCMOS;
                           //double vTh = tech.vth;
                           double vTh = 0.3;
                           double vT = 0.026; // Thermal voltage at room temp (kT/q)

                           // 4. Calculate "On" State Resistance (Linear Region of Tr)
                           if (vGateOn > vTh) {
                               cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
                           } else {
                               cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
                           }

                           // 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
                           double n_factor = 1.2; // Subthreshold swing factor
                           double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));

                           cell.resMemCellOff = cell.readVoltage / iSubOff;

              } else if (cell.mem_rdo == Type::dro) {

		           // 1. In DR, the read voltage is high enough to flip the domain, dumping
                           // the full polarization charge onto the gate. This is represented with a massive capFeOn.

                           double cGateTr = 0.017e-15;

                           // DR Parameters
                           double capFeOn_DR = 80e-15;  // Massive effective switching capacitance
                           double capFeOff = 2e-15;       // Unswitched state (dielectric only)

                           // 2. Calculate Voltage Divider
                           double vGateOn = cell.readVoltage * (capFeOn_DR / (capFeOn_DR + cGateTr));
                           double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));

                           // 3. Extract Transistor parameters for physics equations

                           //Add custom values, tech is not the same
                           //double muCox = tech.mobility * tech.cox;
                           double muCox = 150e-6;
                           double w_l = cell.widthAccessCMOS;
                           //double vTh = tech.vth;
                           double vTh = 0.3;
                           double vT = 0.026; // Thermal voltage at room temp (kT/q)

                           // 4. Calculate "On" State Resistance (Linear Region of Tr)
                           if (vGateOn > vTh) {
                               cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
                           } else {
                               cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
                           }

                           // 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
                           double n_factor = 1.2; // Subthreshold swing factor
                           double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));

                           cell.resMemCellOff = cell.readVoltage / iSubOff;
              }

                double res_equivalent =  (cell.resMemCellOff * cell.resMemCellOn) / (cell.resMemCellOff + cell.resMemCellOn);
		double resTg = CalculateOnResistance(tech.featureSize * 2, NMOS, inputParameter.temperature, tech);

		// 3. Initialize Drivers

        	// WL (Wordline) - Physical Rows
        	wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
                wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

        	// PL (Plateline) - Cols (Physical)
        	plDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numCol)), false, false);
                plSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

        	// BL (Bitline) - Plane (Rows)
                blDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(bitsPerCell)), false, false);
                blSwitchMatrix.Initialize(ROW_MODE, bitsPerCell, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);


        	// Initialize Matrices
                //double resTg = cell.resMemCellAvg / numRow;

                if (cell.writeVoltage > 1.5) {
                        //BL Level Shifter (Row Plane - Horizontal)
                        blLevelShifter.Initialize(bitsPerCell, activityRowWrite, clkFreq);

                        //WL Level Shifter (Row)
                        wlLevelShifter.Initialize(numRow, activityRowWrite, clkFreq);

                        //PL Level Shifter (Column)
                        plLevelShifter.Initialize(numCol, activityColWrite, clkFreq);
                }

		
		if (numColMuxed>1) {
                                //mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, res_equivalent, FPGA);
                                mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);
                                muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
                        }
                        if (param->SARADC) {
                                sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
                        } else {
                                multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
                                multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
                        }

                        currentSenseAmp.Initialize(numCol, false, false, clkFreq, 1);

                        if (numCellPerSynapse > 1) {
                                shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
                        }
                        if (numReadPulse > 1) {
                                shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
                        }

			if (numAdd > 1) {
                                int adderBit = log2(levelOutput) + ceil(log2(numAdd));
                                int numAdder = ceil(numCol/numColMuxed);
                                dff.Initialize(adderBit*numAdder, clkFreq);
                                adder.Initialize(adderBit-1, numAdder, clkFreq);
                        }


		
        	//precharger.Initialize(numCol, resCol, activityColWrite, numReadCellPerOperationNeuro, numWriteCellPerOperationNeuro);
	
	} else if (cell.memCellType == Type::_1T1C) {
        		// Standard 2D 1T1C Array (No 3D Stacking)
        		double cellHeight = cell.heightInFeatureSize;
        		double cellWidth = cell.widthInFeatureSize;

        		if (relaxArrayCellWidth) {
        		    lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;
        		} else {
        		    lengthRow = (double)numCol * cellWidth * tech.featureSize;
        		}
        		if (relaxArrayCellHeight) {
        		    lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
        		} else {
        		    lengthCol = (double)numRow * cellHeight * tech.featureSize;
        		}

        		param->arraywidthunit = cellWidth * tech.featureSize;
        		param->arrayheight = (double)numRow * cellHeight * tech.featureSize;

        		// Parasitics
        		double gateCap = CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech);
        		double drainCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);
        		double sourceCap = CalculateDrainCap(cell.widthAccessCMOS * tech.featureSize, NMOS, cell.widthInFeatureSize * tech.featureSize, tech);

        		// WL (Rows - Gates)
        		capWL = lengthRow * 0.2e-15/1e-6 + (gateCap * numCol);
        		// BL (Cols - Drains/Capacitors)
        		capBL = lengthCol * 0.2e-15/1e-6 + (drainCap * numRow);

        		resRow = lengthRow * unitWireRes;
        		resCol = lengthCol * unitWireRes;

        		param->unitcap = capBL / param->numColSubArray;
        		param->unitres = resCol / param->numRowSubArray;
        		param->columncap = capBL;

        		// Tr Resistance Calculations
        		double cGateTr = gateCap;
        		double capFeOn = 10e-15;  
        		double capFeOff = 2e-15;  
        		double vGateOn = cell.readVoltage * (capFeOn / (capFeOn + cGateTr));
        		double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));
        		
        		double muCox = 150e-6;
        		double w_l = cell.widthAccessCMOS;
        		double vTh = 0.3;
        		double vT = 0.026; 

        		if (vGateOn > vTh) {
        		    cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
        		} else {
        		    cell.resMemCellOn = 1e9; 
        		}
        		double n_factor = 1.2; 
        		double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));
        		cell.resMemCellOff = cell.readVoltage / iSubOff;

        		// Initialize Peripheral Drivers (Reduced count for 1T1C)
        		double resTg = CalculateOnResistance(tech.featureSize * 2, NMOS, inputParameter.temperature, tech);

        		// WL (Row Mode)
        		wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
        		wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

        		// BL (Col Mode - Vertical down the array)
        		blDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numCol)), false, false);
        		blSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

        		if (cell.writeVoltage > 1.5) {
        		    wlLevelShifter.Initialize(numRow, activityRowWrite, clkFreq);
        		    blLevelShifter.Initialize(numCol, activityColWrite, clkFreq);
        		}

        		// Readout Circuits
        		if (numColMuxed > 1) {
        		    mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);
        		    muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
        		}
        		if (param->SARADC) {
        		    sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
        		} else {
        		    multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
        		    multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
        		}
        		currentSenseAmp.Initialize(numCol, false, false, clkFreq, 1);

        		if (numCellPerSynapse > 1) shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
        		if (numReadPulse > 1) shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
        		if (numAdd > 1) {
        		    int adderBit = log2(levelOutput) + ceil(log2(numAdd));
        		    int numAdder = ceil(numCol/numColMuxed);
        		    dff.Initialize(adderBit*numAdder, clkFreq);
        		    adder.Initialize(adderBit-1, numAdder, clkFreq);
        		}
	} else if (cell.memCellType == Type::RRAM ||  cell.memCellType == Type::FeFET || cell.memCellType == Type::Cap) {  //if array is RRAM
		// nvCap added
		double cellHeight = cell.heightInFeatureSize; 
		double cellWidth = cell.widthInFeatureSize;  
		if (cell.accessType == CMOS_access) {  // 1T1R
			if (relaxArrayCellWidth) {
				lengthRow = (double)numCol * MAX(cellWidth, MIN_CELL_WIDTH*2) * tech.featureSize;	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
			} else {
				lengthRow = (double)numCol * cellWidth * tech.featureSize;
			}
			if (relaxArrayCellHeight) {
				lengthCol = (double)numRow * MAX(cellHeight, MIN_CELL_HEIGHT) * tech.featureSize;
			} else {
				lengthCol = (double)numRow * cellHeight * tech.featureSize;
			}
		} else {	// Cross-point, if enter anything else except 'CMOS_access'
			if (relaxArrayCellWidth) {
				lengthRow = (double)numCol * MAX(cellWidth*cell.featureSize, MIN_CELL_WIDTH*2*tech.featureSize);	// Width*2 because generally switch matrix has 2 pass gates per column, even the SL/BL driver has 2 pass gates per column in traditional 1T1R memory
			} else {
				lengthRow = (double)numCol * cellWidth * cell.featureSize;
			}
			if (relaxArrayCellHeight) {
				lengthCol = (double)numRow * MAX(cellHeight*cell.featureSize, MIN_CELL_HEIGHT*tech.featureSize);
			} else {  
				lengthCol = (double)numRow * cellHeight * cell.featureSize;
			}
		}

		// 230920 update
		param->arraywidthunit = cellWidth * cell.featureSize;
		param->arrayheight = (double)numRow * cellHeight * cell.featureSize;

	}      //finish setting array size
	
	// 1.4 update
	capRow1 = lengthRow * 0.2e-15/1e-6;	// BL for 1T1R, WL for Cross-point and SRAM
	capRow2 = lengthRow * 0.2e-15/1e-6;	// WL for 1T1R
	capCol = lengthCol * 0.2e-15/1e-6;

	resRow = lengthRow * param->Metal1_unitwireresis; 
	resCol = lengthCol * param->Metal0_unitwireresis;
	
	
	param->columncap = capCol;
	//start to initializing the subarray modules
	if (cell.memCellType == Type::SRAM) {  //if array is SRAM
		
		//firstly calculate the CMOS resistance and capacitance

		// 1.4 update : modified the code - no folding for SRAM

		// 1.4 update: needs check  - capCol 

		resCellAccess = CalculateOnResistance(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech);

		capCellAccess = CalculateDrainCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, MAX_TRANSISTOR_HEIGHT * tech.featureSize, tech);
		cell.capSRAMCell = capCellAccess + CalculateDrainCap(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, MAX_TRANSISTOR_HEIGHT * tech.featureSize, tech) 
						+ CalculateDrainCap(cell.widthSRAMCellPMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, PMOS, MAX_TRANSISTOR_HEIGHT * tech.featureSize, tech) 
						+ CalculateGateCap(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) + CalculateGateCap(cell.widthSRAMCellPMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech);

		// 1.4 update: for buffer insertion
		double unitcap= capRow1/param->numColSubArray;
		double unitres= resRow/param->numColSubArray;
		param->unitcap = unitcap;
		param->unitres = unitres;	

		if (tech.featureSize <= 14 * 1e-9) capCol += tech.cap_draintotal * cell.widthAccessCMOS * tech.effective_width * numRow;
		else capCol += CalculateDrainCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, MAX_TRANSISTOR_HEIGHT * tech.featureSize, tech) * numRow;	
		param->columncap = capCol;
		if (conventionalSequential) {
		// 1.4 update: consider SRAM parasitic cap
		capRow1 += 2*CalculateGateCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap

			wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
			senseAmp.Initialize(numCol, false, cell.minSenseVoltage, lengthRow/numCol, clkFreq, numReadCellPerOperationNeuro);
			int adderBit = (int)ceil(log2(numRow)) + 1;				
			int numAdder = numCol/numCellPerSynapse; 
			// Anni update: no mux, so adder and shiftaddweight should be for all columns
			dff.Initialize(adderBit*numCol, clkFreq);	
			adder.Initialize(adderBit-1, numCol, clkFreq);
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, numCellPerSynapse);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numAdder, adderBit + numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
			}
			
		} else if (conventionalParallel) {
		// 1.4 update: consider SRAM parasitic cap - only one WL is activated for ADC 
		capRow1 += CalculateGateCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap

			// 1.4 update : buffer insertion
			if (param->buffernumber > 0) {

				sectionres = resRow / (param->buffernumber +1);
				targetdriveres = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize*param->buffersizeratio, NMOS, inputParameter.temperature, tech) ;

				widthInvN  = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, NMOS, inputParameter.temperature, tech) / targetdriveres * tech.featureSize;
				widthInvP = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, PMOS, inputParameter.temperature, tech) / targetdriveres * tech.featureSize ;

				if (tech.featureSize <= 14*1e-9){
					widthInvN = 2* ceil(widthInvN/tech.featureSize) * tech.featureSize;
					widthInvP = 2* ceil(widthInvP/tech.featureSize) * tech.featureSize;
				}				
					

				wlSwitchMatrix.Initialize(ROW_MODE, numRow, sectionres, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			} else {								  
				wlSwitchMatrix.Initialize(ROW_MODE, numRow, resRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			} 
			if (numColMuxed>1) {
				// 1.4 update: store access resistance for multilevelsenseamp
				param->resCellAccess=resCellAccess;
				// 1.4 update: half turn on assume;	Anni update: numRow->numRowParallel
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resCellAccess/(numRowParallel/2), FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}		
			// Anni update: add partial sums	
			int adderBit = log2(levelOutput) + ceil(log2(numAdd));	
			int numAdder = ceil(numCol/numColMuxed);
			if (numAdd > 1) {				
				dff.Initialize(adderBit*numAdder, clkFreq);	
				adder.Initialize(adderBit-1, numAdder, clkFreq);
			}	
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, numCellPerSynapse);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numAdder, adderBit + numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
			}					
		} else if (BNNsequentialMode || XNORsequentialMode) {
			// 1.4 update: consider SRAM parasitic cap
			capRow1 += 2*CalculateGateCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap

			wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
			senseAmp.Initialize(numCol, false, cell.minSenseVoltage, lengthRow/numCol, clkFreq, numReadCellPerOperationNeuro);
			// Anni update
			int adderBit = (int)ceil(log2(numRow)) + 1;	
			int numAdder = numCol;
			dff.Initialize(adderBit*numAdder, clkFreq);	
			adder.Initialize(adderBit-1, numAdder, clkFreq);
		} else if (BNNparallelMode || XNORparallelMode) {
		// 1.4 update: consider SRAM parasitic cap - only one WL is activated for ADC 
		capRow1 += CalculateGateCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap

			wlSwitchMatrix.Initialize(ROW_MODE, numRow, resRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			// Anni update: add mux
			if (numColMuxed>1) {
				// 1.4 update: half turn on assume;	Anni update: numRow->numRowParallel
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resCellAccess/(numRowParallel/2), FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}
			// Anni update: add partial sums
			int adderBit = log2(levelOutput) + ceil(log2(numAdd));	
			int numAdder = ceil(numCol/numColMuxed);
			if (numAdd > 1) {				
				dff.Initialize(adderBit*numAdder, clkFreq);	
				adder.Initialize(adderBit-1, numAdder, clkFreq);
			}
		}
		precharger.Initialize(numCol, resCol, activityColWrite, numReadCellPerOperationNeuro, numWriteCellPerOperationNeuro);
		sramWriteDriver.Initialize(numCol, activityColWrite, numWriteCellPerOperationNeuro);
		
    } else if (cell.memCellType == Type::Cap) { // nvCap added
		cell.resMemCellOn = cell.resistanceOn;        //calculate single memory cell resistance_ON
		cell.resMemCellOff = cell.resistanceOff;      //calculate single memory cell resistance_OFF
		cell.resMemCellAvg = 1/(1/(cell.resistanceOn) * numRowParallel/2.0 + 1/(cell.resistanceOff)* numRowParallel/2.0) * numRowParallel;      //calculate single memory cell resistance_AVG
		if (cell.writeVoltage > 1.5) {
			wllevelshifter.Initialize(numRow, activityRowRead, clkFreq);
			bllevelshifter.Initialize(numRow, activityRowRead, clkFreq);
			sllevelshifter.Initialize(numCol, activityColWrite, clkFreq);
		}
		if (conventionalParallel) { 
			// 1.4 update: needs check enabled rows?;	Anni update: numRow -> numRowParallel
			double resTg = cell.resMemCellAvg / numRowParallel;
			wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg*numRow/numCol, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg * numRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     
			if (numColMuxed>1) {
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}
			// Anni update: add partial sums
			int adderBit = log2(levelOutput) + ceil(log2(numAdd));
			int numAdder = ceil(numCol/numColMuxed);
			if (numAdd > 1) {				
				dff.Initialize(adderBit*numAdder, clkFreq);	
				adder.Initialize(adderBit-1, numAdder, clkFreq);
			}
			// Anni update: shift avgWeightBit for (numCellPerSynapse-1) times, +1 to compensate -1 inside shift-add which designed for general 1-bit cell
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, (numCellPerSynapse-1)*avgWeightBit+1);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numAdder, adderBit + (numCellPerSynapse-1)*avgWeightBit+1, clkFreq, spikingMode, numReadPulse);
			}									
		}	
	


	} else if (cell.memCellType == Type::_2TnC) {


                        // cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
                        // cell.widthAccessCMOS = CalculateOnResistance(tech.featureSize, NMOS, 300, tech) * LINEAR_REGION_RATIO / cell.resCellAccess;   //get access CMOS width
                        // if (cell.widthAccessCMOS > cell.widthInFeatureSize) {   // Place transistor vertically
                        //         printf("Transistor width of 1T1R=%.2fF is larger than the assigned cell width=%.2fF in layout\n", cell.widthAccessCMOS, cell.widthInFeatureSize);
                        //         exit(-1);
                        // }

                        // cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;        //calculate single memory cell resistance_ON
                        // cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF
                        // cell.resMemCellAvg = cell.resCellAccess + cell.resistanceAvg;      //calculate single memory cell resistance_AVG


    			// 2T-nC READ OPERATION PHYSICS (Voltage Coupling & Tr Resistance)
    			//if (cell.memCellType == Type::_2TnC) {

				if (cell.mem_rdo == Type::ndro) {
    			    
    			    		// 1. Calculate the capacitive load at the internal node (Gate of Tr)
    			    		//double cGateTr = CalculateGateCap(cell.widthAccessCMOS * tech.featureSize, tech);
    			    		double cGateTr = 0.017e-15;

    			    		// Define FeCap parameters
    			    		// These represent the switching capacitance differences.
    			    		double capFeOn = 10e-15;  // High switching capacitance (State '0')
    			    		double capFeOff = 2e-15;  // Low switching capacitance (State '1')

    			    		// 2. Voltage Coupling: Calculate V_int (Gate Voltage) for both states
    			    		double vGateOn = cell.readVoltage * (capFeOn / (capFeOn + cGateTr));
    			    		double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));

    			    		// 3. Extract Transistor parameters for physics equations
			    		
			    		//Add custom values, tech is not the same
    			    		//double muCox = tech.mobility * tech.cox;
			    		double muCox = 150e-6;
    			    		double w_l = cell.widthAccessCMOS;
    			    		//double vTh = tech.vth;
    			    		double vTh = 0.3;
    			    		double vT = 0.026; // Thermal voltage at room temp (kT/q)

    			    		// 4. Calculate "On" State Resistance (Linear Region of Tr)
    			    		if (vGateOn > vTh) {
    			    		    cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
    			    		} else {
    			    		    cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
    			    		}

    			    		// 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
    			    		double n_factor = 1.2; // Subthreshold swing factor
    			    		double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));
    			    		
    			    		cell.resMemCellOff = cell.readVoltage / iSubOff;
				
				} else if (cell.mem_rdo == Type::qndro) {
					
					// 1. In QNDR, the read voltage is high enough to flip the domain, dumping 
					// the full polarization charge onto the gate. This is represented with a massive capFeOn.
					
					double cGateTr = 0.017e-15;
					
					// QNDR Parameters
					double capFeOn_QNDR = 40e-15;  // Massive effective switching capacitance
					double capFeOff = 2e-15;       // Unswitched state (dielectric only)
					
					// 2. Calculate Voltage Divider
					double vGateOn = cell.readVoltage * (capFeOn_QNDR / (capFeOn_QNDR + cGateTr));
					double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));

					// 3. Extract Transistor parameters for physics equations

                                        //double muCox = tech.mobility * tech.cox;
                                        double muCox = 150e-6;
                                        double w_l = cell.widthAccessCMOS;
                                        //double vTh = tech.vth;
                                        double vTh = 0.3;
                                        double vT = 0.026; // Thermal voltage at room temp (kT/q)

                                        // 4. Calculate "On" State Resistance (Linear Region of Tr)
                                        if (vGateOn > vTh) {
                                            cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
                                        } else {
                                            cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
                                        }

                                        // 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
                                        double n_factor = 1.2; // Subthreshold swing factor
                                        double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));

                                        cell.resMemCellOff = cell.readVoltage / iSubOff;
				}

				 else if (cell.mem_rdo == Type::dro) {
 
                                          // 1. In DR, the read voltage is high enough to flip the domain, dumping
                                          // the full polarization charge onto the gate. This is represented with a massive capFeOn.
 
                                          double cGateTr = 0.017e-15;
 
                                          // DR Parameters
                                          double capFeOn_DR = 80e-15;  // Massive effective switching capacitance
                                          double capFeOff = 2e-15;       // Unswitched state (dielectric only)
 
                                          // 2. Calculate Voltage Divider
                                          double vGateOn = cell.readVoltage * (capFeOn_DR / (capFeOn_DR + cGateTr));
                                          double vGateOff = cell.readVoltage * (capFeOff / (capFeOff + cGateTr));
 
                                          // 3. Extract Transistor parameters for physics equations
 
                                          //double muCox = tech.mobility * tech.cox;
                                          double muCox = 150e-6;
                                          double w_l = cell.widthAccessCMOS;
                                          //double vTh = tech.vth;
                                          double vTh = 0.3;
                                          double vT = 0.026; // Thermal voltage at room temp (kT/q)
 
                                          // 4. Calculate "On" State Resistance (Linear Region of Tr)
                                          if (vGateOn > vTh) {
                                              cell.resMemCellOn = 1.0 / (muCox * w_l * (vGateOn - vTh));
                                          } else {
					      cell.resMemCellOn = 1e9; // Fallback to High-Z if it doesn't turn on
                                          }
 
                                          // 5. Calculate "Off" State Resistance (Subthreshold Leakage of Tr)
                                          double n_factor = 1.2; // Subthreshold swing factor
                                          double iSubOff = muCox * w_l * (vT * vT) * exp((vGateOff - vTh) / (n_factor * vT));
 
                                          cell.resMemCellOff = cell.readVoltage / iSubOff;
                                  }

				 double res_equivalent =  (cell.resMemCellOff * cell.resMemCellOn)/(cell.resMemCellOff + cell.resMemCellOn);

                        // Initialize Matrices
        		//double resTg = cell.resMemCellAvg / numRow;
			double resTg = CalculateOnResistance(tech.featureSize * 2, NMOS, inputParameter.temperature, tech);

                        if (cell.writeVoltage > 1.5) {
                                //WBL Level Shifter (Row Plane - Horizontal)
                                wblLevelShifter.Initialize(bitsPerCell, activityRowWrite, clkFreq);

                                //SSL Level Shifter (Col)
                                sslLevelShifter.Initialize(numCol, activityColWrite, clkFreq);

                                //WWL Level Shifter (Row)
                                wwlLevelShifter.Initialize(numRow, activityRowWrite, clkFreq);

                                //WPL Level Shifter (Row)
                                wplLevelShifter.Initialize(numRow, activityRowWrite, clkFreq);
                }


                        // double resTg = cell.resMemCellOn / numRow;
                        //double resTg = cell.resMemCellAvg / numRow;

                        //WWL: Write Wordline (ROW MODE)
                        wwlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
                        wwlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

                        //WPL: Write Plateline (ROW MODE)
                        wplDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
                        wplSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);


                        //RBL: Read Bitline (ROW MODE)
                        rblDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
                        rblSwitchMatrix.Initialize(ROW_MODE, numRow, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

                        //WBL: Write Bitline (ROW MODE)
                        wblDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);
                        wblPlaneDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(bitsPerCell)), false, false);

                        wblSwitchMatrix.Initialize(ROW_MODE, numRow * bitsPerCell, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);


                        //SSL: Source Select Line (COL MODE)
                        sslDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numCol)), false, false);
                        sslSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);

                        //RSL: Read Select/Source Line (COL MODE)
                        rslDecoder.Initialize(REGULAR_COL, (int)ceil(log2(numCol)), false, false);
                        rslSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);



                        if (numColMuxed>1) {
                                // mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, res_equivalent, FPGA);
                                mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);
                                muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
                        }
                        if (param->SARADC) {
                                sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
                        } else {
                                multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
                                multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
                        }

			currentSenseAmp.Initialize(numCol, false, false, clkFreq, 1);

                        if (numCellPerSynapse > 1) {
                                shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
                        }
                        if (numReadPulse > 1) {
                                shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
                        }

			if (numAdd > 1) {
    			        int adderBit = log2(levelOutput) + ceil(log2(numAdd));
    			        int numAdder = ceil(numCol/numColMuxed);
    			        dff.Initialize(adderBit*numAdder, clkFreq);
    			        adder.Initialize(adderBit-1, numAdder, clkFreq);
    			}




	} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
		if (cell.accessType == CMOS_access) {	// 1T1R

			cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
			cell.widthAccessCMOS = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, NMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO / cell.resCellAccess;   //get access CMOS width
			double widthAccessInFeatureSize = cell.widthAccessCMOS;
			if (tech.featureSize <= 14 * 1e-9){			
				widthAccessInFeatureSize = ((cell.widthAccessCMOS-1) * tech.PitchFin + tech.widthFin) / tech.featureSize;  //convert #fin to F
			}
			if (widthAccessInFeatureSize > cell.widthInFeatureSize) {	// Place transistor vertically
				printf("Transistor width of 1T1R=%.2fF is larger than the assigned cell width=%.2fF in layout\n", cell.widthAccessCMOS, cell.widthInFeatureSize);
				exit(-1);
			}
			cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;        //calculate single memory cell resistance_ON
			cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF

			// 1.4 update; Anni update: numRow->numRowParallel
			cell.resMemCellAvg = 1/(1/(cell.resistanceOn + cell.resCellAccess ) * numRowParallel/2.0 + 1/(cell.resistanceOff + cell.resCellAccess )* numRowParallel/2.0) * numRowParallel;      //calculate single memory cell resistance_AVG

			// 1.4 update: needs check  - capCol : cell.widthInFeatureSize / 
			capRow2 += CalculateGateCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) * numCol;          //sum up all the gate cap of access CMOS, as the row cap
			// 1.4 update
			double unitcap= capRow2/param->numColSubArray;
			double unitres= resRow/param->numColSubArray;
			param->unitcap = unitcap;
			param->unitres = unitres;			
			
			capCol += CalculateDrainCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, MAX_TRANSISTOR_HEIGHT * tech.featureSize, tech) * numRow;	// If capCol is found to be too large, increase cell.widthInFeatureSize to relax the limit
			param->columncap = capCol;
		} else {	// Cross-point
			cell.resMemCellOn = cell.resistanceOn;
			cell.resMemCellOff = cell.resistanceOff;
			cell.resMemCellOnAtHalfVw = cell.resistanceOn;
			cell.resMemCellOffAtHalfVw = cell.resistanceOff;
			cell.resMemCellOnAtVw = cell.resistanceOn;
			cell.resMemCellOffAtVw = cell.resistanceOff;

			// 1.4 update; Anni update: numRow->numRowParallel
			cell.resMemCellAvg = 1/(1/(cell.resistanceOn + cell.resCellAccess ) * numRowParallel/2.0 + 1/(cell.resistanceOff + cell.resCellAccess )* numRowParallel/2.0) * numRowParallel;      //calculate single memory cell resistance_AVG

			cell.resMemCellAvgAtHalfVw = cell.resistanceAvg;
			cell.resMemCellAvgAtVw = cell.resistanceAvg;
		}
		
		if (cell.writeVoltage > 1.5) {
			wllevelshifter.Initialize(numRow, activityRowRead, clkFreq);
			bllevelshifter.Initialize(numRow, activityRowRead, clkFreq);
			sllevelshifter.Initialize(numCol, activityColWrite, clkFreq);
		}
		
		if (conventionalSequential) {  
			double capBL = lengthCol * 0.2e-15/1e-6;
			int numAdder = (int)ceil(numCol/numColMuxed);   // numCol is divisible by numCellPerSynapse
			int numInput = numAdder;        //XXX input number of MUX, 
			double resTg = cell.resMemCellOn;     //transmission gate resistance
			int adderBit = (int)ceil(log2(numRow)) + avgWeightBit;  
						
			wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);          
			if (cell.accessType == CMOS_access) {
				wlNewDecoderDriver.Initialize(numRow);          
			} else {
				wlDecoderDriver.Initialize(ROW_MODE, numRow, numCol);
			}						
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     //SL use switch matrix
			if (numColMuxed>1) {
				mux.Initialize(numInput, numColMuxed, resTg, FPGA);     
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, pow(2, avgWeightBit), clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, pow(2, avgWeightBit), clkFreq, numReadCellPerOperationNeuro, false, currentMode);
				if (avgWeightBit > 1) {
					multilevelSAEncoder.Initialize(pow(2, avgWeightBit), numCol/numColMuxed);
				}
			}
			// Anni update
			dff.Initialize(adderBit*numAdder, clkFreq); 
			adder.Initialize(adderBit-1, numAdder, clkFreq);
			// Anni update: shift avgWeightBit for (numCellPerSynapse-1) times, +1 to compensate -1 inside shift-add which designed for general 1-bit cell
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, (numCellPerSynapse-1)*avgWeightBit+1);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numAdder, adderBit + (numCellPerSynapse-1)*avgWeightBit+1, clkFreq, spikingMode, numReadPulse);
			}
			
		} else if (conventionalParallel) { 
			// 1.4 update: needs check enabled rows?;	Anni update: numRow -> numRowParallel
			double resTg = cell.resMemCellAvg / numRowParallel;
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.Initialize(numRow, activityRowRead, clkFreq);
				// 1.4 update: buffer insertion
				if (param->buffernumber>0) {
					sectionres = resRow / (param->buffernumber +1);
					targetdriveres = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize*param->buffersizeratio, NMOS, inputParameter.temperature, tech) ;

					widthInvN  = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, NMOS, inputParameter.temperature, tech) / targetdriveres * tech.featureSize;
					widthInvP = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, PMOS, inputParameter.temperature, tech) / targetdriveres * tech.featureSize ;

					if (tech.featureSize <= 14*1e-9){
						widthInvN = 2* ceil(widthInvN/tech.featureSize) * tech.featureSize;
						widthInvP = 2* ceil(widthInvP/tech.featureSize) * tech.featureSize;
					}															
				}
			} else {
				wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg*numRow/numCol, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			}
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg * numRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     
			if (numColMuxed>1) {
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}
			// Anni update: add partial sums
			int adderBit = log2(levelOutput) + ceil(log2(numAdd));
			int numAdder = ceil(numCol/numColMuxed);
			if (numAdd > 1) {				
				dff.Initialize(adderBit*numAdder, clkFreq);	
				adder.Initialize(adderBit-1, numAdder, clkFreq);
			}
			// Anni update: shift avgWeightBit for (numCellPerSynapse-1) times, +1 to compensate -1 inside shift-add which designed for general 1-bit cell
			// needs check
			/*
			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(numAdder, adderBit, clkFreq, spikingMode, (numCellPerSynapse-1)*avgWeightBit+1);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(numAdder, adderBit + (numCellPerSynapse-1)*avgWeightBit+1, clkFreq, spikingMode, numReadPulse);
			}	
			*/					

			if (numCellPerSynapse > 1) {
				shiftAddWeight.Initialize(ceil(numCol/numColMuxed), log2(levelOutput), clkFreq, spikingMode, numCellPerSynapse);
			}
			if (numReadPulse > 1) {
				shiftAddInput.Initialize(ceil(numCol/numColMuxed), log2(levelOutput)+numCellPerSynapse, clkFreq, spikingMode, numReadPulse);
			}		

		} else if (BNNsequentialMode || XNORsequentialMode) {       
			double resTg = cell.resMemCellOn;
			int numAdder = (int)ceil(numCol/numColMuxed);  
			int numInput = numAdder;        
			int adderBit = (int)ceil(log2(numRow)) + 1; 
			
			wlDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numRow)), false, false);           
			if (cell.accessType == CMOS_access) {
				wlNewDecoderDriver.Initialize(numRow);          
			} else {
				wlDecoderDriver.Initialize(ROW_MODE, numRow, numCol);
			}
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     //SL use switch matrix
			if (numColMuxed>1) {
				mux.Initialize(numInput, numColMuxed, resTg, FPGA);      
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed)), true, false);
			}

			// 1.4 update 230615
			multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, false, currentMode);
			multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);

			// rowCurrentSenseAmp.Initialize(numCol/numColMuxed, true, false, clkFreq, numReadCellPerOperationNeuro);
			
			// Anni update
			dff.Initialize(adderBit*numAdder, clkFreq); 
			adder.Initialize(adderBit-1, numAdder, clkFreq);
		} else if (BNNparallelMode || XNORparallelMode) {      
			// 1.4 update: needs check enabled rows?;	Anni update: numRow->numRowParallel
			double resTg = cell.resMemCellAvg / numRowParallel;
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.Initialize(numRow, activityRowRead, clkFreq);         
			} else {
				wlSwitchMatrix.Initialize(ROW_MODE, numRow, resTg*numRow/numCol, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);
			}
			slSwitchMatrix.Initialize(COL_MODE, numCol, resTg*numRow, true, false, activityRowRead, activityColWrite, numWriteCellPerOperationMemory, numWriteCellPerOperationNeuro, 1, clkFreq);     
			if (numColMuxed>1) {
				mux.Initialize(ceil(numCol/numColMuxed), numColMuxed, resTg, FPGA);       
				muxDecoder.Initialize(REGULAR_ROW, (int)ceil(log2(numColMuxed/2)), true, true);    
			}
			if (param->SARADC) {
				sarADC.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro);
			} else {
				multilevelSenseAmp.Initialize(numCol/numColMuxed, levelOutput, clkFreq, numReadCellPerOperationNeuro, true, currentMode);
				multilevelSAEncoder.Initialize(levelOutput, numCol/numColMuxed);
			}
			// Anni update: add partial sums
			int adderBit = log2(levelOutput) + ceil(log2(numAdd));	
			int numAdder = ceil(numCol/numColMuxed);
			if (numAdd > 1) {				
				dff.Initialize(adderBit*numAdder, clkFreq);	
				adder.Initialize(adderBit-1, numAdder, clkFreq);
			}
		}
	} 
	initialized = true;  //finish initialization

	if (cell.memCellType == Type::_2TnC) {
	    double capFeChk   = EPS0_VACUUM * param->epsFE * param->cellAreaFE / param->tFE;
	    double capGateChk = capFeChk / param->capGateRatio;
	    double capNodeChk = param->bitsPerCell * capFeChk + capGateChk;
	    double dVsenseChk = (param->twoPr * param->cellAreaFE) / capNodeChk;
	    double nMaxChk    = param->twoPr * param->tFE
	                      / (EPS0_VACUUM * param->epsFE * param->senseAmpResolution);
	    cout << "[2T-nC sense] n=" << param->bitsPerCell
	         << "  C_node=" << capNodeChk*1e15 << " fF"
	         << "  dV=" << dVsenseChk*1e3 << " mV"
	         << "  (amp resolves " << param->senseAmpResolution*1e3
	         << " mV, n_max=" << nMaxChk << ")" << endl;
	    if (dVsenseChk < param->senseAmpResolution)
	        cout << " Differential signal below the amplifier's resolution."
	             << " Raise 2Pr or tFE, or use a lower-offset sense amp." << endl;
	}
}



void SubArray::CalculateArea() {  //calculate layout area for total design
	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;  //ensure initialization first
	} else {  //if initialized, start to do calculation
		area = 0;
		usedArea = 0;
		if (cell.memCellType == Type::SRAM) {       
			// Array only
			heightArray = lengthCol;
			widthArray = lengthRow;
			areaArray = heightArray * widthArray;
			
			//precharger and writeDriver are always needed for all different designs
			precharger.CalculateArea(NULL, widthArray, NONE);
			sramWriteDriver.CalculateArea(NULL, widthArray, NONE);
			
			if (conventionalSequential) {
				wlDecoder.CalculateArea(heightArray, NULL, NONE);  
				senseAmp.CalculateArea(NULL, widthArray, MAGIC);
				adder.CalculateArea(NULL, widthArray, NONE);
				dff.CalculateArea(NULL, widthArray, NONE);
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: DFF*2, for adder and shift-add weight pipeline
				height = precharger.height + sramWriteDriver.height + heightArray + senseAmp.height + adder.height + dff.height*2 + shiftAddInput.height + shiftAddWeight.height;
				width = wlDecoder.width + widthArray;
				area = height * width;
				usedArea = areaArray + wlDecoder.area + precharger.area + sramWriteDriver.area + senseAmp.area + adder.area + dff.area*2 + shiftAddInput.area + shiftAddWeight.area;
				emptyArea = area - usedArea;

				areaADC = senseAmp.area + precharger.area;
				areaAccum = adder.area + dff.area*2 + shiftAddInput.area + shiftAddWeight.area;
				areaOther = wlDecoder.area + sramWriteDriver.area;
			} else if (conventionalParallel) { 
				wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				if (numColMuxed>1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: add partial sums, height & usedArea & areaAccum
				if (numAdd > 1) {
					adder.CalculateArea(NULL, widthArray, NONE);
					dff.CalculateArea(NULL, widthArray, NONE);
				}
				
				// 1.4 update : repeater implementation
				// buffer area
				if (param->buffernumber>0) {					
					CalculateGateArea(INV, 1, widthInvN , widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
					CalculateGateCapacitance(INV, 1, widthInvN , widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &drivecapin, &drivecapout);					
				} else {
					wInv = 0;
					hInv = 0;
				}
				// 1.4 update : repeater implementation
				double bufferarea= hInv * wInv * param->buffernumber * 2 * param->numRowSubArray;
				height = precharger.height + sramWriteDriver.height + heightArray + multilevelSenseAmp.height + multilevelSAEncoder.height + \
						shiftAddInput.height + shiftAddWeight.height + ((numAdd > 1)==true? (adder.height+dff.height):0) + ((numColMuxed > 1)==true? (mux.height):0)+sarADC.height;
				width = MAX(wlSwitchMatrix.width, ((numColMuxed > 1)==true? (muxDecoder.width):0)) + widthArray + bufferarea/lengthCol; // added for buffer area;
				area = height * width;
				usedArea = areaArray + wlSwitchMatrix.area + precharger.area + sramWriteDriver.area + multilevelSenseAmp.area + multilevelSAEncoder.area + \
						shiftAddInput.area + shiftAddWeight.area + ((numAdd > 1)==true? (adder.area+dff.area):0) + ((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0)+sarADC.area + bufferarea;
				emptyArea = area - usedArea;

				areaADC = multilevelSenseAmp.area + precharger.area + multilevelSAEncoder.area + sarADC.area;
				areaAccum = shiftAddInput.area + shiftAddWeight.area + ((numAdd > 1)==true? (adder.area+dff.area):0);
				areaOther = wlSwitchMatrix.area + sramWriteDriver.area + ((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + bufferarea;
			} else if (BNNsequentialMode || XNORsequentialMode) {
				wlDecoder.CalculateArea(heightArray, NULL, NONE);  
				senseAmp.CalculateArea(NULL, widthArray, MAGIC);
				adder.CalculateArea(NULL, widthArray, NONE);
				dff.CalculateArea(NULL, widthArray, NONE);
				height = precharger.height + sramWriteDriver.height + heightArray + senseAmp.height + adder.height + dff.height;
				width = wlDecoder.width + widthArray;
				area = height * width;
				usedArea = areaArray + wlDecoder.area + precharger.area + sramWriteDriver.area + senseAmp.area + adder.area + dff.area;
				emptyArea = area - usedArea;
			} else if (BNNparallelMode || XNORparallelMode) {
				wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				// Anni update: add mux
				if (numColMuxed>1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: add partial sums, height & usedArea
				if (numAdd > 1) {
					adder.CalculateArea(NULL, widthArray, NONE);
					dff.CalculateArea(NULL, widthArray, NONE);
				}
				height = precharger.height + sramWriteDriver.height + heightArray + multilevelSenseAmp.height + multilevelSAEncoder.height + sarADC.height + \
						((numColMuxed > 1)==true? (mux.height):0) + ((numAdd > 1)==true? (adder.height+dff.height):0);
				width = MAX(wlSwitchMatrix.width, ((numColMuxed > 1)==true? (muxDecoder.width):0)) + widthArray;
				area = height * width;
				usedArea = areaArray + wlSwitchMatrix.area + precharger.area + sramWriteDriver.area + multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area + \
							((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + ((numAdd > 1)==true? (adder.area+dff.area):0);
				emptyArea = area - usedArea;
			}
	    } 
		else if (cell.memCellType == Type::Cap) { //nvCap added
			// v1.5 update: floorplan with less area
			heightArray = lengthCol;
			widthArray = lengthRow;
			areaArray = heightArray * widthArray;
			
			if (conventionalParallel) { 
				
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSenseAmp.area /= 10;	// area for VSA only excluding voltage divider
					multilevelSenseAmp.area += 0.5*numCol/numColMuxed*5*656e-15; // area for opamp
					multilevelSenseAmp.height = multilevelSenseAmp.area/multilevelSenseAmp.width; // height for opamp
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: add partial sums, height & usedArea & areaAccum
				if (numAdd > 1) {
					adder.CalculateArea(NULL, widthArray, NONE);
					dff.CalculateArea(NULL, widthArray, NONE);
				}
				
				height = ((cell.writeVoltage > 1.5)==true? (sllevelshifter.height):0) + slSwitchMatrix.height + heightArray + ((numColMuxed > 1)==true? (mux.height):0) + \
						multilevelSenseAmp.height + multilevelSAEncoder.height + shiftAddWeight.height + shiftAddInput.height + ((numAdd > 1)==true? (adder.height+dff.height):0) + sarADC.height;
				
				wlSwitchMatrix.CalculateArea(height-muxDecoder.height, NULL, NONE);
				wllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				bllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);

				width = MAX( ((cell.writeVoltage > 1.5)==true? (wllevelshifter.width + bllevelshifter.width):0) + wlSwitchMatrix.width, ((numColMuxed > 1)==true? (muxDecoder.width):0)) + widthArray;
				usedArea = areaArray + ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlSwitchMatrix.area + slSwitchMatrix.area + 
							((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + multilevelSenseAmp.area  + multilevelSAEncoder.area + shiftAddWeight.area + shiftAddInput.area + ((numAdd > 1)==true? (adder.area+dff.area):0) + sarADC.area;
				
				areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
				areaAccum = shiftAddWeight.area + shiftAddInput.area + ((numAdd > 1)==true? (adder.area+dff.area):0);
				areaOther = ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlSwitchMatrix.area + slSwitchMatrix.area + ((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0);
				
				area = height * width;				
				emptyArea = area - usedArea;

			}
			else{
				cout << "[Subarray] Error: Conventional parallel mode is required for nvCap!" << endl;
				exit(-1);
			}
		


		} else if (cell.memCellType == Type::_2TnC) {
                double bitsPerCell = param->bitsPerCell;
		//Calculate Core Array Area
                heightArray = lengthCol;
                widthArray = lengthRow;
                areaArray = heightArray * widthArray;

                //Peripheral Circuits (Assumed similar to FeFET/DRAM)
               //  wblDecoder.CalculateArea(heightArray, NULL, NONE);
               //  wblSwitchMatrix.CalculateArea(heightArray, NULL, NONE);

                //SSL (Source Select Line)
                sslDecoder.CalculateArea(NULL, widthArray, NONE);
                sslSwitchMatrix.CalculateArea(NULL, widthArray, NONE);

                //RSL (Read Select Line)
                rslDecoder.CalculateArea(NULL, widthArray, NONE);
                rslSwitchMatrix.CalculateArea(NULL, widthArray, NONE);


               //  //WWL (Write Wordline)
               //  wwlDecoder.CalculateArea(heightArray, NULL, NONE);
               //  wwlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);

               //  //WPL (Write Plateline)
               //  wplDecoder.CalculateArea(heightArray, NULL, NONE);
               //  wplSwitchMatrix.CalculateArea(heightArray, NULL, NONE);

               //  //RBL (Read Bitline - Standard)
               //  rblDecoder.CalculateArea(heightArray, NULL, NONE);
               //  rblSwitchMatrix.CalculateArea(heightArray, NULL, NONE);


                //Mux & ADC
                if (numColMuxed > 1) {
                    mux.CalculateArea(NULL, widthArray, NONE);
                    muxDecoder.CalculateArea(NULL, NULL, NONE);
                    double minMuxHeight = MAX(muxDecoder.height, mux.height);
                    mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
                }

		// if (cell.writeVoltage > 1.5) {
            	//     sslLevelShifter.CalculateArea(NULL, widthArray, NONE);
            	// }

                if (param->SARADC) {
                        // SAR ADC calculates its natural area without constraints
                        sarADC.CalculateUnitArea();
                        sarADC.CalculateArea(NULL, widthArray, NONE);
                } else {
                        // Flash ADC (MLSA)
                        multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
                        multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
                }

		//currentSenseAmp.CalculateArea(widthArray);
		currentSenseAmp.CalculateUnitArea();
		currentSenseAmp.CalculateArea(widthArray/numCol);


                //Shift-Add Logic (Accumulators)
                if (numReadPulse > 1) {
                        shiftAddInput.CalculateArea(NULL, widthArray, NONE);
                }
                if (numCellPerSynapse > 1) {
                        shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
                }
		if (numAdd > 1) {
    		        adder.CalculateArea(NULL, widthArray, NONE);
    		        dff.CalculateArea(NULL, widthArray, NONE);
    		}


                //Total Height = Array + Column Drivers + Readout Circuits
		//Add level shifters here
                height = heightArray
                                 + rslSwitchMatrix.height + rslDecoder.height
                                 + sslSwitchMatrix.height + sslDecoder.height
                                 // + ((cell.writeVoltage > 1.5) ? sslLevelShifter.height : 0)
				 + ((numColMuxed > 1) ? mux.height : 0)
                                 + (param->SARADC ? sarADC.height : (multilevelSenseAmp.height + multilevelSAEncoder.height))
				 // + currentSenseAmp.height
                                 + ((numReadPulse > 1) ? shiftAddInput.height : 0)
				 + ((numCellPerSynapse > 1) ? shiftAddWeight.height : 0)
				 + ((numAdd > 1) ? (adder.height + dff.height) : 0);

		//double rowDriverHeight = height - ((numColMuxed > 1) ? muxDecoder.height : 0);
		
				double rowDriverHeight = height - ((numColMuxed > 1) ? muxDecoder.height : 0);

		{
			double minBlk = MAX(ceil(log2((double)numRow)), 1.0)
			              * tech.featureSize * MAX_TRANSISTOR_HEIGHT;
			if (rowDriverHeight < minBlk) rowDriverHeight = minBlk;
		}


		if (cell.writeVoltage > 1.5) {
                    // Row-oriented Level Shifters (Constrained by heightArray)
                    wwlLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);
                    wplLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);
                    wblLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);

                    // Column-oriented Level Shifter (Constrained by widthArray)
                    //sslLevelShifter.CalculateArea(NULL, widthArray, NONE);
                }

                wblDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
		wblPlaneDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
		wblSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

		//WWL (Write Wordline)
                wwlDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
                wwlSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

                //WPL (Write Plateline)
                wplDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
                wplSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

                //RBL (Read Bitline - Standard)
                rblDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
                rblSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

		double totalRowDriverWidth = ((cell.writeVoltage > 1.5) ? (wwlLevelShifter.width + wplLevelShifter.width + wblLevelShifter.width) : 0)
		            + wwlSwitchMatrix.width + wwlDecoder.width
                            + wplSwitchMatrix.width + wplDecoder.width
                            + wblSwitchMatrix.width + wblDecoder.width + wblPlaneDecoder.width
                            + rblSwitchMatrix.width + rblDecoder.width;


		//Total Width = Array + Row Drivers
                width = widthArray + 
                		MAX(totalRowDriverWidth, ((numColMuxed > 1) ? muxDecoder.width : 0));


                //Calculate Final Area
                area = height * width;

                //Calculate Used Area (Active silicon only, excluding empty space/white space)
		usedArea = areaArray
                       + ((cell.writeVoltage > 1.5) ? (wwlLevelShifter.area + wplLevelShifter.area + wblLevelShifter.area) : 0)
		       + wwlSwitchMatrix.area + wwlDecoder.area
                       + wplSwitchMatrix.area + wplDecoder.area
                       + rblSwitchMatrix.area + rblDecoder.area
                       + wblSwitchMatrix.area + wblDecoder.area + wblPlaneDecoder.area
                       + sslSwitchMatrix.area + sslDecoder.area
                       + rslSwitchMatrix.area + rslDecoder.area
                       + ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
                       + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area))
		       + currentSenseAmp.area
                       + ((numReadPulse > 1) ? shiftAddInput.area : 0)
		       + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
	               + ((numAdd > 1) ? (adder.area + dff.area) : 0); 

                if (area < usedArea) { cout << "wtf THIS IS NOT CORRECT" << endl ;}
		emptyArea = area - usedArea;

		areaADC = (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area));
		areaADC = areaADC + currentSenseAmp.area;
                areaAccum = ((numReadPulse > 1) ? shiftAddInput.area : 0) + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0) + ((numAdd > 1) ? (adder.area + dff.area) : 0);
                areaOther = usedArea - areaArray - areaADC - areaAccum;

//		{
//			int    nLayers3D     = (int)param->bitsPerCell;   // stacked WBL planes
//			double staircaseStep = 400e-9;                    // landing pitch per plane [ASSUMPTION]
//			double depthPerSide  = ceil(nLayers3D / 2.0) * staircaseStep;
//
//			areaStaircase = 2.0 * depthPerSide * heightArray;
//			areaCore      = areaArray + areaStaircase;
//
//			// the two strips are dead silicon appended to the width edges
//			width    += 2.0 * depthPerSide;
//			area      = height * width;
		{
			int    nLayers3D     = (int)param->bitsPerCell;
			double staircaseStep = param->staircaseStep;      // was 400e-9 local 
			double padPitch      = param->staircasePadPitch;  // along the strip  

			int    perSide = (int)ceil(nLayers3D / 2.0);
			int    nx      = (int)ceil(sqrt((double)perSide));   // steps in depth 
			int    ny      = (int)ceil((double)perSide / nx);    // rows of pads   

			double depthPerSide = nx * staircaseStep;
			//double stripSpan    = MIN(heightArray, ny * padPitch);
			double stripSpan = MIN(MIN(heightArray, widthArray), ny * padPitch);

			areaStaircase = 2.0 * depthPerSide * stripSpan;
			areaCore      = areaArray + areaStaircase;

			width    += 2.0 * depthPerSide;
			area      = height * width;
			emptyArea = area - usedArea;	
		}

		// 3D floorplan: staircase + architecture
	//	{
	//	int nLayers3D = (int)param->bitsPerCell;   // stacked WBL planes

 	//       StaircaseResult stair = CalculateDoubleSidedStaircase(nLayers3D, heightArray);
 	//       areaStaircase  = stair.staircaseArea;
	//       areaCore = areaArray + areaStaircase;
 	//       areaViaRouting = stair.viaRoutingArea;
 	//       areaBondPad    = 0.0;


 	//       // READOUT GROUP -> placed directly UNDER THE ARRAY.
 	//       // Everything the RSL columns land on, plus the SSL data drivers
 	//       // (column-pitched, write path) and the accumulation chain.
 	//       double areaReadoutGroup =
 	//       	  currentSenseAmp.area
 	//       	+ ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
 	//       	+ sslSwitchMatrix.area + sslDecoder.area
 	//       	+ ((numReadPulse > 1)      ? shiftAddInput.area  : 0)
 	//       	+ ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
 	//       	+ ((numAdd > 1)            ? (adder.area + dff.area) : 0);
 	//       areaReadoutGroup += (param->SARADC ? sarADC.area
 	//                                          : (multilevelSenseAmp.area + multilevelSAEncoder.area));
 	//       areaReadoutGroup += rslSwitchMatrix.area + rslDecoder.area;   // to be removed

 	//       double areaDriverGroup =
 	//       	  wwlSwitchMatrix.area + wwlDecoder.area
 	//       	+ wplSwitchMatrix.area + wplDecoder.area
 	//       	+ rblSwitchMatrix.area + rblDecoder.area
 	//       	+ wblSwitchMatrix.area + wblDecoder.area
 	//       	+ ((cell.writeVoltage > 1.5) ? (wwlLevelShifter.area + wplLevelShifter.area + wblLevelShifter.area) : 0);

 	//       usedArea  = areaArray + areaReadoutGroup + areaDriverGroup;
 	//       areaADC   = currentSenseAmp.area
 	//                 + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area));
 	//       areaAccum = ((numReadPulse > 1) ? shiftAddInput.area : 0)
 	//                 + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
 	//                 + ((numAdd > 1) ? (adder.area + dff.area) : 0);
 	//       areaOther = usedArea - areaArray - areaADC - areaAccum;

 	//       if (param->integrationMode == CNA) {
 	//       	width        += 2.0 * stair.depthPerSide;
 	//       	area          = height * width;
 	//       	emptyArea     = area - usedArea;
 	//       	areaMemoryDie = area;
 	//       	areaLogicDie  = 0.0;
 	//       	chipFootprint = area;

 	//       } else {
 	//       	StackedDieResult die = CalculateStackedDie(
 	//       		widthArray, heightArray, stair,
 	//       		areaReadoutGroup, areaDriverGroup, areaArray,
 	//       		param->integrationMode, param->numBondPadPerSubarray);

 	//       	areaMemoryDie = die.memoryDieArea;
 	//       	areaLogicDie  = die.logicDieArea;
 	//       	areaBondPad   = die.bondPadArea;
 	//       	areaViaRouting += die.tavArea;
 	//       	chipFootprint = die.chipFootprint;

 	//       	area      = die.chipFootprint;    // top-down outline (dies overlap in Z)
 	//       	usedArea  = die.activeSilicon;
 	//       	emptyArea = die.deadArea;

 	//       	height = heightArray;
	//	width  = area / height;

	//       }
	       
//}	

		
		// cout << "\n=================== 2TnC Area Breakdown (um^2) ===================" << endl;
		// cout << "Array Core Area:           " << areaArray * 1e12 << endl;
		// cout << "WWL Decoder Area:          " << wwlDecoder.area * 1e12 << endl;
		// cout << "WWL Switch Matrix Area:    " << wwlSwitchMatrix.area * 1e12 << endl;
		// cout << "WPL Decoder Area:          " << wplDecoder.area * 1e12 << endl;
		// cout << "WPL Switch Matrix Area:    " << wplSwitchMatrix.area * 1e12 << endl;
		// cout << "WBL Decoder Area:          " << wblDecoder.area * 1e12 << endl;
		// cout << "WBL Switch Matrix Area:    " << wblSwitchMatrix.area * 1e12 << endl;
		// cout << "RSL Decoder Area:          " << rslDecoder.area * 1e12 << endl;
		// cout << "RSL Switch Matrix Area:    " << rslSwitchMatrix.area * 1e12 << endl;
		// cout << "RBL Decoder Area:          " << rblDecoder.area * 1e12 << endl;
		// cout << "RBL Switch Matrix Area:    " << rblSwitchMatrix.area * 1e12 << endl;
		// if (cell.writeVoltage > 1.5) {
		//     cout << "WWL Level Shifter Area:    " << wwlLevelShifter.area * 1e12 << endl;
		//     cout << "WPL Level Shifter Area:    " << wplLevelShifter.area * 1e12 << endl;
		//     cout << "WBL Level Shifter Area:    " << wblLevelShifter.area * 1e12 << endl;
		// }
		// if (numColMuxed > 1) {
		//     cout << "MUX Area:                  " << mux.area * 1e12 << endl;
		//     cout << "MUX Decoder Area:          " << muxDecoder.area * 1e12 << endl;
		// }
		// if (param->SARADC) {
		//     cout << "SAR ADC Area:              " << sarADC.area * 1e12 << endl;
		// } else {
		//     cout << "Multilevel SA Area:        " << multilevelSenseAmp.area * 1e12 << endl;
		//     cout << "Multilevel SA Encoder Area:" << multilevelSAEncoder.area * 1e12 << endl;
		// }
		// cout << "Current Sense Amp Area:    " << currentSenseAmp.area * 1e12 << endl;
		// if (numReadPulse > 1) cout << "ShiftAdd Input Area:       " << shiftAddInput.area * 1e12 << endl;
		// if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Area:      " << shiftAddWeight.area * 1e12 << endl;
		// if (numAdd > 1) {
		//     cout << "Adder Area:                " << adder.area * 1e12 << endl;
		//     cout << "DFF Area:                  " << dff.area * 1e12 << endl;
		// }
		// cout << "------------------------------------------------------------------" << endl;
		// cout << "Used Area (Active Silicon):" << usedArea * 1e12 << endl;
		// cout << "Empty Area (White Space):  " << emptyArea * 1e12 << endl;
		// cout << "TOTAL SUBARRAY AREA:       " << area * 1e12 << endl;
		// cout << "==================================================================\n" << endl;

		} else if (cell.memCellType == Type::_1TnC) {
        		double bitsPerCell = param->bitsPerCell;

			heightArray = lengthCol;
        		widthArray = lengthRow;
        		areaArray = heightArray * widthArray;
			// cout << "Raw Memory Cell Area (No Peripherals): " << areaArray * 1e12 << " um^2" << endl;

			plSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
                        plDecoder.CalculateArea(NULL, widthArray, NONE);

        		//Mux & ADC
               		if (numColMuxed > 1) {
               		     mux.CalculateArea(NULL, widthArray, NONE);
               		     muxDecoder.CalculateArea(NULL, NULL, NONE);
               		     double minMuxHeight = MAX(muxDecoder.height, mux.height);
               		     mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
               		}

			if (cell.writeVoltage > 1.5) {
				plLevelShifter.CalculateArea(NULL, widthArray, NONE);
			}

			if (param->SARADC) {
                        // SAR ADC calculates its natural area without constraints
                        sarADC.CalculateUnitArea();
                        sarADC.CalculateArea(NULL, widthArray, NONE);
                	} else {
                	        // Flash ADC (MLSA)
                	        multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
                	        multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
                	}

                	currentSenseAmp.CalculateUnitArea();
			//currentSenseAmp.CalculateArea(widthArray);
			currentSenseAmp.CalculateArea(widthArray/numCol);


                	//Shift-Add Logic (Accumulators)
                	if (numReadPulse > 1) {
                	        shiftAddInput.CalculateArea(NULL, widthArray, NONE);
                	}
                	if (numCellPerSynapse > 1) {
                	        shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
                	}
                	if (numAdd > 1) {
                	        adder.CalculateArea(NULL, widthArray, NONE);
                	        dff.CalculateArea(NULL, widthArray, NONE);
                	}

			//Total Height = Array + Column Drivers + Readout Circuits
                	//Add level shifters here
                	// double maxColPeripheralWidth = 0;
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, plSwitchMatrix.width);
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, plDecoder.width);
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, ((numColMuxed > 1) ? mux.width : 0));
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, (param->SARADC ? sarADC.width : (multilevelSenseAmp.width + multilevelSAEncoder.width)));
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, currentSenseAmp.width);
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, ((numReadPulse > 1) ? shiftAddInput.width : 0));
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, ((numCellPerSynapse > 1) ? shiftAddWeight.width : 0));
        		// maxColPeripheralWidth = MAX(maxColPeripheralWidth, ((numAdd > 1) ? MAX(adder.width, dff.width) : 0));
        		// if (cell.writeVoltage > 1.5) maxColPeripheralWidth = MAX(maxColPeripheralWidth, plLevelShifter.width);

			// double effectiveArrayWidth = MAX(widthArray, maxColPeripheralWidth);

			// // Calculate the TRUE height of the row drivers
        		// if (cell.writeVoltage > 1.5) {
        		//     wlLevelShifter.CalculateArea(heightArray, NULL, NONE);
        		//     blLevelShifter.CalculateArea(heightArray, NULL, NONE);
        		// }
        		// blDecoder.CalculateArea(heightArray, NULL, NONE);
        		// blSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
        		// wlDecoder.CalculateArea(heightArray, NULL, NONE);
        		// wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);

        		// double maxRowPeripheralHeight = 0;
        		// maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, blDecoder.height);
        		// maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, blSwitchMatrix.height);
        		// maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, wlDecoder.height);
        		// maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, wlSwitchMatrix.height);
        		// if (cell.writeVoltage > 1.5) {
        		//     maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, wlLevelShifter.height);
        		//     maxRowPeripheralHeight = MAX(maxRowPeripheralHeight, blLevelShifter.height);
        		// }

			// double effectiveArrayHeight = MAX(heightArray, maxRowPeripheralHeight);

			// // Total Height = Effective Array Height + Column Drivers
        		// height = effectiveArrayHeight
        		//          + ((cell.writeVoltage > 1.5) ? plLevelShifter.height : 0)
        		//          + plSwitchMatrix.height + plDecoder.height
        		//          + ((numColMuxed > 1) ? mux.height : 0)
        		//          + (param->SARADC ? sarADC.height : (multilevelSenseAmp.height + multilevelSAEncoder.height))
        		//          + currentSenseAmp.height
        		//          + ((numReadPulse > 1) ? shiftAddInput.height : 0)
        		//          + ((numCellPerSynapse > 1) ? shiftAddWeight.height : 0)
        		//          + ((numAdd > 1) ? (adder.height + dff.height) : 0);

        		// // Total Width = Effective Array Width + Row Drivers
        		// double totalRowDriverWidth = ((cell.writeVoltage > 1.5) ? (wlLevelShifter.width + blLevelShifter.width) : 0)
        		//                              + wlSwitchMatrix.width + wlDecoder.width
        		//                              + blSwitchMatrix.width + blDecoder.width;

        		// width = effectiveArrayWidth + MAX(totalRowDriverWidth, ((numColMuxed > 1) ? muxDecoder.width : 0));
			
			 height = heightArray
                	                  + ((cell.writeVoltage > 1.5) ? plLevelShifter.height : 0)
                	                  + plSwitchMatrix.height + plDecoder.height
                	                  + ((numColMuxed > 1) ? mux.height : 0)
                	                  + (param->SARADC ? sarADC.height : (multilevelSenseAmp.height + multilevelSAEncoder.height))
                	                  + ((numReadPulse > 1) ? shiftAddInput.height : 0)
			 		  + ((numCellPerSynapse > 1) ? shiftAddWeight.height : 0)
                	                  + ((numAdd > 1) ? (adder.height + dff.height) : 0);

			 double rowDriverHeight = height - ((numColMuxed > 1) ? muxDecoder.height : 0);

			 if (cell.writeVoltage > 1.5) {
                	     // Row-oriented Level Shifters (Constrained by heightArray)
                	     wlLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);
                	     blLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);
                	 }

                	 blDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
                	 blSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

                	 //WL (Wordline)
                	 wlDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
                	 wlSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

			 double totalRowDriverWidth = ((cell.writeVoltage > 1.5) ? (wlLevelShifter.width + blLevelShifter.width) : 0)
               		 		             + wlSwitchMatrix.width + wlDecoder.width
               		 		             + blSwitchMatrix.width + blDecoder.width;


               		 //Total Width = Array + Row Drivers
               		 width = widthArray + MAX(totalRowDriverWidth, ((numColMuxed > 1) ? muxDecoder.width : 0));


			area = height * width;

			//Calculate Used Area (Active silicon only, excluding empty space/white space)
               		 // usedArea = areaArray
               		 //        + ((cell.writeVoltage > 1.5) ? (wlLevelShifter.area + plLevelShifter.area + blLevelShifter.area) : 0)
               		 //        + wlSwitchMatrix.area + wlDecoder.area
               		 //        + plSwitchMatrix.area + plDecoder.area
               		 //        + blSwitchMatrix.area + blDecoder.area
               		 //        + ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
               		 //        + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area))
               		 //        + currentSenseAmp.area
			 //        + ((numReadPulse > 1) ? shiftAddInput.area : 0)
	                 //        + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
               		 //        + ((numAdd > 1) ? (adder.area + dff.area) : 0);

			 usedArea = areaArray
               		        + ((cell.writeVoltage > 1.5) ? (wlLevelShifter.area + plLevelShifter.area + blLevelShifter.area) : 0)
               		        + wlSwitchMatrix.area + wlDecoder.area
               		        + plSwitchMatrix.area + plDecoder.area
               		        + blSwitchMatrix.area + blDecoder.area
               		        + ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
               		        + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area))
				+ ((numReadPulse > 1) ? shiftAddInput.area : 0)
	                        + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
               		        + ((numAdd > 1) ? (adder.area + dff.area) : 0);

                     // Safeguard against Mux Decoder corner overlap
                     if (area < usedArea) { cout << "wtf THIS IS NOT CORRECT" << endl ;}
               		 emptyArea = area - usedArea;

               		 //emptyArea = area - usedArea;

               		 areaADC = (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area));
			 areaAccum = ((numReadPulse > 1) ? shiftAddInput.area : 0) + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0) + ((numAdd > 1) ? (adder.area + dff.area) : 0);
               		 areaOther = usedArea - areaArray - areaADC - areaAccum;

			 {
				int    nLayers3D     = (int)param->bitsPerCell;
				double staircaseStep = 400e-9;
				double depthPerSide  = ceil(nLayers3D / 2.0) * staircaseStep;

				areaStaircase = 2.0 * depthPerSide * heightArray;
				areaCore      = areaArray + areaStaircase;

				width    += 2.0 * depthPerSide;
				area      = height * width;
				emptyArea = area - usedArea;
			 }
			 
			 //  // 3D floorplan: staircase + architecture 
			//  {
			// 	int nLayers3D = (int)param->bitsPerCell;   // stacked plate (BL) planes

			// 	StaircaseResult stair = CalculateDoubleSidedStaircase(nLayers3D, heightArray);
			// 	areaStaircase  = stair.staircaseArea;
			// 	areaCore = areaArray + areaStaircase;
			// 	areaViaRouting = stair.viaRoutingArea;
			// 	areaBondPad    = 0.0;

			// 	double areaReadoutGroup =
			// 		  currentSenseAmp.area
			// 		+ ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
			// 		+ plSwitchMatrix.area + plDecoder.area
			// 		+ ((numReadPulse > 1)      ? shiftAddInput.area  : 0)
			// 		+ ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
			// 		+ ((numAdd > 1)            ? (adder.area + dff.area) : 0)
			// 		+ (param->SARADC ? sarADC.area
			// 		                 : (multilevelSenseAmp.area + multilevelSAEncoder.area));

			// 	double areaDriverGroup =
			// 		  wlSwitchMatrix.area + wlDecoder.area
			// 		+ blSwitchMatrix.area + blDecoder.area
			// 		+ ((cell.writeVoltage > 1.5) ? (wlLevelShifter.area + blLevelShifter.area + plLevelShifter.area) : 0);

			// 	usedArea  = areaArray + areaReadoutGroup + areaDriverGroup;
			// 	areaADC   = currentSenseAmp.area
			// 	          + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area));
			// 	areaAccum = ((numReadPulse > 1) ? shiftAddInput.area : 0)
			// 	          + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
			// 	          + ((numAdd > 1) ? (adder.area + dff.area) : 0);
			// 	areaOther = usedArea - areaArray - areaADC - areaAccum;

			// 	if (param->integrationMode == CNA) {
			// 		width        += 2.0 * stair.depthPerSide;
			// 		area          = height * width;
			// 		emptyArea     = area - usedArea;
			// 		areaMemoryDie = area;
			// 		areaLogicDie  = 0.0;
			// 		chipFootprint = area;
			// 	} else {
			// 		StackedDieResult die = CalculateStackedDie(
			// 			widthArray, heightArray, stair,
			// 			areaReadoutGroup, areaDriverGroup, areaArray,
			// 			param->integrationMode, param->numBondPadPerSubarray);

			// 		areaMemoryDie  = die.memoryDieArea;
			// 		areaLogicDie   = die.logicDieArea;
			// 		areaBondPad    = die.bondPadArea;
			// 		areaViaRouting += die.tavArea;
			// 		chipFootprint  = die.chipFootprint;

			// 		area      = die.chipFootprint;
			// 		usedArea  = die.activeSilicon;
			// 		emptyArea = die.deadArea;
			// 		height    = heightArray;
			// 		width     = area / height;
			// 	}

//}

			 //  cout << "\n=================== 1TnC Area Breakdown (um^2) ===================" << endl;
		 	 // cout << "Array Core Area:           " << areaArray * 1e12 << endl;
		 	 // cout << "WL Decoder Area:           " << wlDecoder.area * 1e12 << endl;
		 	 // cout << "WL Switch Matrix Area:     " << wlSwitchMatrix.area * 1e12 << endl;
		 	 // cout << "BL Decoder Area:           " << blDecoder.area * 1e12 << endl;
		 	 // cout << "BL Switch Matrix Area:     " << blSwitchMatrix.area * 1e12 << endl;
		 	 // cout << "PL Decoder Area:           " << plDecoder.area * 1e12 << endl;
		 	 // cout << "PL Switch Matrix Area:     " << plSwitchMatrix.area * 1e12 << endl;
		 	 // if (cell.writeVoltage > 1.5) {
		 	 //     cout << "WL Level Shifter Area:     " << wlLevelShifter.area * 1e12 << endl;
		 	 //     cout << "BL Level Shifter Area:     " << blLevelShifter.area * 1e12 << endl;
		 	 //     cout << "PL Level Shifter Area:     " << plLevelShifter.area * 1e12 << endl;
		 	 // }
		 	 // if (numColMuxed > 1) {
		 	 //     cout << "MUX Area:                  " << mux.area * 1e12 << endl;
		 	 //     cout << "MUX Decoder Area:          " << muxDecoder.area * 1e12 << endl;
		 	 // }
		 	 // if (param->SARADC) {
		 	 //     cout << "SAR ADC Area:              " << sarADC.area * 1e12 << endl;
		 	 // } else {
		 	 //     cout << "Multilevel SA Area:        " << multilevelSenseAmp.area * 1e12 << endl;
		 	 //     cout << "Multilevel SA Encoder Area:" << multilevelSAEncoder.area * 1e12 << endl;
		 	 // }
		 	 // if (numReadPulse > 1) cout << "ShiftAdd Input Area:       " << shiftAddInput.area * 1e12 << endl;
		 	 // if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Area:      " << shiftAddWeight.area * 1e12 << endl;
		 	 // if (numAdd > 1) {
		 	 //     cout << "Adder Area:                " << adder.area * 1e12 << endl;
		 	 //     cout << "DFF Area:                  " << dff.area * 1e12 << endl;
		 	 // }
		 	 // cout << "------------------------------------------------------------------" << endl;
		 	 // cout << "Used Area (Active Silicon):" << usedArea * 1e12 << endl;
		 	 // cout << "Empty Area (White Space):  " << emptyArea * 1e12 << endl;
		 	 // cout << "TOTAL SUBARRAY AREA:       " << area * 1e12 << endl;
		 	 // cout << "==================================================================\n" << endl;


		} else if (cell.memCellType == Type::_1T1C) {
        		heightArray = lengthCol;
        		widthArray = lengthRow;
        		areaArray = heightArray * widthArray;

        		// Column drivers
        		blSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
        		blDecoder.CalculateArea(NULL, widthArray, NONE);
        		if (cell.writeVoltage > 1.5) blLevelShifter.CalculateArea(NULL, widthArray, NONE);

        		if (numColMuxed > 1) {
        		    mux.CalculateArea(NULL, widthArray, NONE);
        		    muxDecoder.CalculateArea(NULL, NULL, NONE);
        		    mux.CalculateArea(MAX(muxDecoder.height, mux.height), widthArray, OVERRIDE);
        		}

        		if (param->SARADC) {
        		    sarADC.CalculateUnitArea();
        		    sarADC.CalculateArea(NULL, widthArray, NONE);
        		} else {
        		    multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
        		    multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
        		}

        		currentSenseAmp.CalculateUnitArea();
        		currentSenseAmp.CalculateArea(widthArray/numCol);

        		if (numReadPulse > 1) shiftAddInput.CalculateArea(NULL, widthArray, NONE);
        		if (numCellPerSynapse > 1) shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
        		if (numAdd > 1) {
        		    adder.CalculateArea(NULL, widthArray, NONE);
        		    dff.CalculateArea(NULL, widthArray, NONE);
        		}

        		// Total Height = Array + Col Drivers + Sense/Math
        		height = heightArray
        		         + ((cell.writeVoltage > 1.5) ? blLevelShifter.height : 0)
        		         + blSwitchMatrix.height + blDecoder.height
        		         + ((numColMuxed > 1) ? mux.height : 0)
        		         + (param->SARADC ? sarADC.height : (multilevelSenseAmp.height + multilevelSAEncoder.height))
        		         + ((numReadPulse > 1) ? shiftAddInput.height : 0)
        		         + ((numCellPerSynapse > 1) ? shiftAddWeight.height : 0)
        		         + ((numAdd > 1) ? (adder.height + dff.height) : 0);

        		double rowDriverHeight = height - ((numColMuxed > 1) ? muxDecoder.height : 0);

        		// Row drivers
        		if (cell.writeVoltage > 1.5) wlLevelShifter.CalculateArea(rowDriverHeight, NULL, NONE);
        		wlDecoder.CalculateArea(rowDriverHeight, NULL, NONE);
        		wlSwitchMatrix.CalculateArea(rowDriverHeight, NULL, NONE);

        		double totalRowDriverWidth = ((cell.writeVoltage > 1.5) ? wlLevelShifter.width : 0)
        		                             + wlSwitchMatrix.width + wlDecoder.width;

        		width = widthArray + MAX(totalRowDriverWidth, ((numColMuxed > 1) ? muxDecoder.width : 0));
        		area = height * width;

        		usedArea = areaArray
        		         + ((cell.writeVoltage > 1.5) ? (wlLevelShifter.area + blLevelShifter.area) : 0)
        		         + wlSwitchMatrix.area + wlDecoder.area
        		         + blSwitchMatrix.area + blDecoder.area
        		         + ((numColMuxed > 1) ? (mux.area + muxDecoder.area) : 0)
        		         + (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area))
        		         + ((numReadPulse > 1) ? shiftAddInput.area : 0)
        		         + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0)
        		         + ((numAdd > 1) ? (adder.area + dff.area) : 0);

        		emptyArea = area - usedArea;

        		areaADC = (param->SARADC ? sarADC.area : (multilevelSenseAmp.area + multilevelSAEncoder.area));
        		areaAccum = ((numReadPulse > 1) ? shiftAddInput.area : 0) + ((numCellPerSynapse > 1) ? shiftAddWeight.area : 0) + ((numAdd > 1) ? (adder.area + dff.area) : 0);
        		areaOther = usedArea - areaArray - areaADC - areaAccum;

		} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			// Array only
			heightArray = lengthCol;
			widthArray = lengthRow;
			areaArray = heightArray * widthArray;
			
			// Level shifter for write
			if (cell.writeVoltage > 1.5) {
				wllevelshifter.CalculateArea(heightArray, NULL, NONE);
				bllevelshifter.CalculateArea(heightArray, NULL, NONE);
				sllevelshifter.CalculateArea(NULL, widthArray, NONE);				
			}
			
			if (conventionalSequential) { 			
				//20241031 update move to the bottom	
				// wlDecoder.CalculateArea(heightArray, NULL, NONE);
				// if (cell.accessType == CMOS_access) {
				// 	wlNewDecoderDriver.CalculateArea(heightArray, NULL, NONE);
				// } else {
				// 	wlDecoderDriver.CalculateArea(heightArray, NULL, NONE);
				// }				
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					if (avgWeightBit > 1) {
						multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
					}
				}
				
				dff.CalculateArea(NULL, widthArray, NONE);
				adder.CalculateArea(NULL, widthArray, NONE);
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
				}
				//20241031 update area floor plannning with less empty area
				height = ((cell.writeVoltage > 1.5)==true? (sllevelshifter.height):0) + slSwitchMatrix.height + heightArray + ((numColMuxed > 1)==true? (mux.height):0) + \
						multilevelSenseAmp.height + multilevelSAEncoder.height + adder.height + dff.height + shiftAddInput.height + shiftAddWeight.height + sarADC.height;
				if (cell.accessType == CMOS_access) {
					wlNewDecoderDriver.CalculateArea(height-muxDecoder.height, NULL, NONE);
				} 
				else {
					wlDecoderDriver.CalculateArea(height-muxDecoder.height, NULL, NONE);
				}
				wllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				bllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				wlDecoder.CalculateArea(height-muxDecoder.height, NULL, NONE);
				width = MAX( ((cell.writeVoltage > 1.5)==true? (wllevelshifter.width + bllevelshifter.width):0) + wlDecoder.width + wlNewDecoderDriver.width + wlDecoderDriver.width, ((numColMuxed > 1)==true? (muxDecoder.width):0) ) + widthArray;
				area = height * width;
				usedArea = areaArray + ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlDecoder.area + wlDecoderDriver.area + wlNewDecoderDriver.area + slSwitchMatrix.area + 
							((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + multilevelSenseAmp.area + multilevelSAEncoder.area + adder.area + dff.area + shiftAddInput.area + shiftAddWeight.area + sarADC.area;
				emptyArea = area - usedArea;
				
				areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
				areaAccum = adder.area + dff.area + shiftAddInput.area + shiftAddWeight.area;
				areaOther = ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlDecoder.area + wlNewDecoderDriver.area + wlDecoderDriver.area + slSwitchMatrix.area + ((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0);
			
			// cout << "\n================= RRAM (Crossbar - Sequential) Area Breakdown (um^2) =================" << endl;
			// 	 cout << "Array Core Area:           " << areaArray * 1e12 << endl;
			// 	 cout << "WL Decoder Area:           " << wlDecoder.area * 1e12 << endl;
			// 	 
			// 	 // Differentiate based on the access type (1T1R vs Cross-point)
			// 	 if (cell.accessType == CMOS_access) {
			// 	     cout << "WL Decoder Driver Area:    " << wlNewDecoderDriver.area * 1e12 << endl;
			// 	 } else {
			// 	     cout << "WL Decoder Driver Area:    " << wlDecoderDriver.area * 1e12 << endl;
			// 	 }
			// 	 
			// 	 cout << "SL Switch Matrix Area:     " << slSwitchMatrix.area * 1e12 << endl;
			// 	 
			// 	 if (cell.writeVoltage > 1.5) {
			// 	     cout << "WL Level Shifter Area:     " << wllevelshifter.area * 1e12 << endl;
			// 	     cout << "BL Level Shifter Area:     " << bllevelshifter.area * 1e12 << endl;
			// 	     cout << "SL Level Shifter Area:     " << sllevelshifter.area * 1e12 << endl;
			// 	 }
			// 	 if (numColMuxed > 1) {
			// 	     cout << "MUX Area:                  " << mux.area * 1e12 << endl;
			// 	     cout << "MUX Decoder Area:          " << muxDecoder.area * 1e12 << endl;
			// 	 }
			// 	 if (param->SARADC) {
			// 	     cout << "SAR ADC Area:              " << sarADC.area * 1e12 << endl;
			// 	 } else {
			// 	     cout << "Multilevel SA Area:        " << multilevelSenseAmp.area * 1e12 << endl;
			// 	     if (avgWeightBit > 1) {
			// 	         cout << "Multilevel SA Encoder Area:" << multilevelSAEncoder.area * 1e12 << endl;
			// 	     }
			// 	 }
			// 	 
			// 	 if (numReadPulse > 1) cout << "ShiftAdd Input Area:       " << shiftAddInput.area * 1e12 << endl;
			// 	 if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Area:      " << shiftAddWeight.area * 1e12 << endl;
			// 	 cout << "Adder Area:                " << adder.area * 1e12 << endl;
			// 	 cout << "DFF Area:                  " << dff.area * 1e12 << endl;
			// 	 cout << "-------------------------------------------------------------------------" << endl;
			// 	 cout << "Used Area (Active Silicon):" << usedArea * 1e12 << endl;
			// 	 cout << "Empty Area (White Space):  " << emptyArea * 1e12 << endl;
			// 	 cout << "TOTAL SUBARRAY AREA:       " << area * 1e12 << endl;
			// 	 cout << "=========================================================================\n" << endl;
			} else if (conventionalParallel) {
				//20241031 update switch to the bottom  
				// if (cell.accessType == CMOS_access) {
				// 	wlNewSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				// } else {
				// 	wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				// }
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculateArea(NULL, widthArray, NONE);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: add partial sums, height & usedArea & areaAccum
				if (numAdd > 1) {
					adder.CalculateArea(NULL, widthArray, NONE);
					dff.CalculateArea(NULL, widthArray, NONE);
				}
				
				// 1.4 update : buffer area
				// buffer area
				if (param->buffernumber>0) {
					CalculateGateArea(INV, 1, widthInvN , widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &hInv, &wInv);
					CalculateGateCapacitance(INV, 1, widthInvN , widthInvP, tech.featureSize * MAX_TRANSISTOR_HEIGHT, tech, &drivecapin, &drivecapout);					
				} else {
					wInv = 0;
					hInv = 0;
				}
				double bufferarea= hInv * wInv * param->buffernumber * 2 * param->numRowSubArray;
				// 20241031 update the floor plan calculation is changed with less empty area.
				height = ((cell.writeVoltage > 1.5)==true? (sllevelshifter.height):0) + slSwitchMatrix.height + heightArray + ((numColMuxed > 1)==true? (mux.height):0) + \
						multilevelSenseAmp.height + multilevelSAEncoder.height + shiftAddWeight.height + shiftAddInput.height + ((numAdd > 1)==true? (adder.height+dff.height):0) + sarADC.height;
				if (cell.accessType == CMOS_access) {
					wlNewSwitchMatrix.CalculateArea(height-muxDecoder.height, NULL, NONE);
				} 
				else {
					wlSwitchMatrix.CalculateArea(height-muxDecoder.height, NULL, NONE);
				}
				wllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				bllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				width = MAX( ((cell.writeVoltage > 1.5)==true? (wllevelshifter.width + bllevelshifter.width):0) + wlNewSwitchMatrix.width, ((numColMuxed > 1)==true? (muxDecoder.width):0)) + widthArray + bufferarea/lengthCol; // added;
				usedArea = areaArray + ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlNewSwitchMatrix.area + slSwitchMatrix.area + 
						((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + multilevelSenseAmp.area  + multilevelSAEncoder.area + shiftAddWeight.area + shiftAddInput.area + ((numAdd > 1)==true? (adder.area+dff.area):0) + sarADC.area + bufferarea;
				
				areaADC = multilevelSenseAmp.area + multilevelSAEncoder.area + sarADC.area;
				areaAccum = shiftAddWeight.area + shiftAddInput.area + ((numAdd > 1)==true? (adder.area+dff.area):0);
				areaOther = ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlNewSwitchMatrix.area + slSwitchMatrix.area + ((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + bufferarea;
				
				area = height * width;				
				emptyArea = area - usedArea;

				// cout << "\n================= RRAM (Crossbar) Area Breakdown (um^2) =================" << endl;
				//  cout << "Array Core Area:           " << areaArray * 1e12 << endl;

                 		// // In Parallel Mode, we use Switch Matrices instead of Decoders
				//  if (cell.accessType == CMOS_access) {
				//      cout << "WL Switch Matrix Area:     " << wlNewSwitchMatrix.area * 1e12 << endl;
				//  } else {
				//      cout << "WL Switch Matrix Area:     " << wlSwitchMatrix.area * 1e12 << endl;
				//  }

                 		// cout << "SL Switch Matrix Area:     " << slSwitchMatrix.area * 1e12 << endl;

                 		// if (cell.writeVoltage > 1.5) {
				//      cout << "WL Level Shifter Area:     " << wllevelshifter.area * 1e12 << endl;
				//      cout << "BL Level Shifter Area:     " << bllevelshifter.area * 1e12 << endl;
				//      cout << "SL Level Shifter Area:     " << sllevelshifter.area * 1e12 << endl;
				//  }
				//  if (numColMuxed > 1) {
				//      cout << "MUX Area:                  " << mux.area * 1e12 << endl;
				//      cout << "MUX Decoder Area:          " << muxDecoder.area * 1e12 << endl;
				//  }
				//  if (param->SARADC) {
				//      cout << "SAR ADC Area:              " << sarADC.area * 1e12 << endl;
				//  } else {
				//      cout << "Multilevel SA Area:        " << multilevelSenseAmp.area * 1e12 << endl;
				//      cout << "Multilevel SA Encoder Area:" << multilevelSAEncoder.area * 1e12 << endl;
				//  }
				//  if (numReadPulse > 1) cout << "ShiftAdd Input Area:       " << shiftAddInput.area * 1e12 << endl;
				//  if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Area:      " << shiftAddWeight.area * 1e12 << endl;
				//  if (numAdd > 1) {
				//      cout << "Adder Area:                " << adder.area * 1e12 << endl;
				//      cout << "DFF Area:                  " << dff.area * 1e12 << endl;
				//  }
				//  cout << "-------------------------------------------------------------------------" << endl;
				//  cout << "Used Area (Active Silicon):" << usedArea * 1e12 << endl;
				//  cout << "Empty Area (White Space):  " << emptyArea * 1e12 << endl;
				//  cout << "TOTAL SUBARRAY AREA:       " << area * 1e12 << endl;
				//  cout << "=========================================================================\n" << endl;


			} else if (BNNsequentialMode || XNORsequentialMode) {    
				wlDecoder.CalculateArea(heightArray, NULL, NONE);
				if (cell.accessType == CMOS_access) {
					wlNewDecoderDriver.CalculateArea(heightArray, NULL, NONE);
				} else {
					wlDecoderDriver.CalculateArea(heightArray, NULL, NONE);
				}
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}

				// 1.4 update 230615

				// rowCurrentSenseAmp.CalculateUnitArea();
				// rowCurrentSenseAmp.CalculateArea(widthArray);
				multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);	
				multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);

				dff.CalculateArea(NULL, widthArray, NONE);
				adder.CalculateArea(NULL, widthArray, NONE);
				
				//20241031 update subarray floor planning with less empty area
				height = ((cell.writeVoltage > 1.5)==true? (sllevelshifter.height):0) + slSwitchMatrix.height + heightArray + ((numColMuxed > 1)==true? (mux.height):0) + multilevelSAEncoder.height + multilevelSenseAmp.height + adder.height + dff.height;
				if (cell.accessType == CMOS_access) {
					wlNewDecoderDriver.CalculateArea(height-muxDecoder.height, NULL, NONE);
				} 
				else {
					wlDecoderDriver.CalculateArea(height-muxDecoder.height, NULL, NONE);
				}
				wllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				bllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				wlDecoder.CalculateArea(height-muxDecoder.height, NULL, NONE);
				width = MAX( ((cell.writeVoltage > 1.5)==true? (wllevelshifter.width + bllevelshifter.width):0) + wlDecoder.width + wlNewDecoderDriver.width + wlDecoderDriver.width, ((numColMuxed > 1)==true? (muxDecoder.width):0)) + widthArray;
				area = height * width;
				// 1.4 update 230615s
				usedArea = areaArray + ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlDecoder.area + wlDecoderDriver.area + wlNewDecoderDriver.area + slSwitchMatrix.area + 
							((numColMuxed > 1)==true? (mux.area + muxDecoder.area):0) + multilevelSenseAmp.area + multilevelSAEncoder.area  + adder.area + dff.area;
				emptyArea = area - usedArea;
			} else if (BNNparallelMode || XNORparallelMode) {      
				if (cell.accessType == CMOS_access) {
					wlNewSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				} else {
					wlSwitchMatrix.CalculateArea(heightArray, NULL, NONE);
				}
				slSwitchMatrix.CalculateArea(NULL, widthArray, NONE);
				if (numColMuxed > 1) {
					mux.CalculateArea(NULL, widthArray, NONE);
					muxDecoder.CalculateArea(NULL, NULL, NONE);
					double minMuxHeight = MAX(muxDecoder.height, mux.height);
					mux.CalculateArea(minMuxHeight, widthArray, OVERRIDE);
				}
				if (param->SARADC) {
					sarADC.CalculateUnitArea();
					sarADC.CalculateArea(NULL, widthArray, NONE);
				} else {
					multilevelSenseAmp.CalculateArea(NULL, widthArray, NONE);
					multilevelSAEncoder.CalculateArea(NULL, widthArray, NONE);
				}
				// Anni update: add partial sums, height & usedArea
				if (numAdd > 1) {
					adder.CalculateArea(NULL, widthArray, NONE);
					dff.CalculateArea(NULL, widthArray, NONE);
				}
				height = ((cell.writeVoltage > 1.5)==true? (sllevelshifter.height):0) + slSwitchMatrix.height + heightArray + mux.height + multilevelSenseAmp.height + multilevelSAEncoder.height + sarADC.height + ((numAdd > 1)==true? (adder.height+dff.height):0);
				// 20241031 update area floor planning with less empty
				if (cell.accessType == CMOS_access) {
					wlNewSwitchMatrix.CalculateArea(height-muxDecoder.height, NULL, NONE);
				} 
				else {
					wlSwitchMatrix.CalculateArea(height-muxDecoder.height, NULL, NONE);
				}
				wllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				bllevelshifter.CalculateArea(height-muxDecoder.height, NULL, NONE);
				width = MAX( ((cell.writeVoltage > 1.5)==true? (wllevelshifter.width + bllevelshifter.width):0) + wlNewSwitchMatrix.width + wlSwitchMatrix.width, muxDecoder.width) + widthArray;
				area = height * width;
				usedArea = areaArray + ((cell.writeVoltage > 1.5)==true? (wllevelshifter.area + bllevelshifter.area + sllevelshifter.area):0) + wlSwitchMatrix.area + wlNewSwitchMatrix.area + slSwitchMatrix.area + 
							mux.area + multilevelSenseAmp.area + muxDecoder.area + multilevelSAEncoder.area + sarADC.area + ((numAdd > 1)==true? (adder.area+dff.area):0);
				emptyArea = area - usedArea;
			}
		} 
	}
}

void SubArray::CalculateLatency(double columnRes, const vector<double> &columnResistance, bool CalculateclkFreq, bool writeBack) {   //calculate latency for different mode 

	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		
		readLatency = 0;
		readLatencyADC = 0;
		readLatencyAccum = 0;
		readLatencyOther = 0;
		writeLatency = 0;

		if (cell.memCellType == Type::SRAM) {
			if (conventionalSequential) {
				int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {

					// 1.4 update: new arguments for rowdecoder
					wlDecoder.CalculateLatency(1e20, capRow1, NULL, resRow, numCol, 1, numRow*activityRowWrite);				
					precharger.CalculateLatency(1e20, capCol, 1, numWriteOperationPerRow*numRow*activityRowWrite);					
					senseAmp.CalculateLatency(1);

					// Read
					// 1.4 update: SRAM column cap update (calibration with IMEC data)
					// 1.4 update: SRAM row delay update
					// needs check

					double resPullDown = CalculateOnResistance(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech);
					double BLCap_perCell = capCol / numRow + capCellAccess; // Anni update
					double BLRes_perCell = resCol / numRow;
					double Elmore_BL = (resCellAccess + resPullDown) * BLCap_perCell * numRow   + BLCap_perCell * BLRes_perCell * numRow  * ( numRow +1 )  /2;

					colDelay = Elmore_BL * log(tech.vdd / (tech.vdd - cell.minSenseVoltage / 2));  

					if (CalculateclkFreq) {
						readLatency += wlDecoder.readLatency;
						readLatency += precharger.readLatency;
						readLatency += colDelay;
						readLatency += senseAmp.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				} 
				if (!CalculateclkFreq) {
					
					// Anni update: hide readLatencyAccum by pipeline
					adder.CalculateLatency(1e20, dff.capTgDrain, 1);	// numRead = numReadOperationPerRow*numRow*activityRowRead
					dff.CalculateLatency(1e20, 1); // numRead = numReadOperationPerRow*numRow*activityRowRead
					if (numCellPerSynapse > 1) {
						shiftAddWeight.CalculateLatency(1);	// numCellPerSynapse-1
					}
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(1);					
					}					
					if (param->synchronous) {
						readLatencyADC = numReadOperationPerRow*numRow*activityRowRead;
						// Anni update: hide readLatencyAccum by pipeline
						// adder is pipelined with ADC
						readLatencyAccum += numReadOperationPerRow*(numRow*activityRowRead-1) * (ceil(adder.readLatency*clkFreq)-1);		
						// shiftAddWeight and shiftAddInput are pipelined with ADC+Adder (no mux, only once computation for whole array)
						readLatencyAccum += MAX(ceil((shiftAddWeight.adder.readLatency*(numCellPerSynapse-1)+shiftAddInput.adder.readLatency)*clkFreq)-(readLatencyADC+readLatencyAccum), 0);	
					} else {
						readLatencyADC = (precharger.readLatency + colDelay + senseAmp.readLatency) * numReadOperationPerRow*numRow*activityRowRead * (validated==true? param->beta : 1);		
						readLatencyOther = wlDecoder.readLatency * numRow*activityRowRead * (validated==true? param->beta : 1);
						// Anni update
						readLatencyAccum = MAX(adder.readLatency*(numReadOperationPerRow*numRow*activityRowRead-1) + shiftAddWeight.adder.readLatency*(numCellPerSynapse-1) + \
											shiftAddInput.adder.readLatency - readLatencyADC - readLatencyOther, 0);
					}
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}	
					// // Write (assume the average delay of pullup and pulldown inverter in SRAM cell)
					// double resPull;
					// resPull = (CalculateOnResistance(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, inputParameter.temperature, tech) + CalculateOnResistance(cell.widthSRAMCellPMOS * tech.featureSize, PMOS, inputParameter.temperature, tech)) / 2;    // take average
					// tau = resPull * cell.capSRAMCell;
					// gm = (CalculateTransconductance(cell.widthSRAMCellNMOS * tech.featureSize, NMOS, tech) + CalculateTransconductance(cell.widthSRAMCellPMOS * tech.featureSize, PMOS, tech)) / 2;   // take average
					// beta = 1 / (resPull * gm);
					// sramWriteDriver.CalculateLatency(1e20, capCol, resCol, numWriteOperationPerRow*numRow*activityRowWrite);
					// // writeLatency += horowitz(tau, beta, 1e20, NULL) * numWriteOperationPerRow * numRow * activityRowWrite;
					// // writeLatency += wlDecoder.writeLatency;
					// // writeLatency += precharger.writeLatency;
					// // writeLatency += sramWriteDriver.writeLatency;
			} else if (conventionalParallel) {
				int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {
					double bufferlatency = 0;

					// 1.4 update : buffer for the latency
					int iterbuffnum = param->buffernumber-1;
					if (param->buffernumber>0) {
						double buffnum = param->buffernumber; // # of buffers
						double sectionnum = param->numColSubArray/(param->buffernumber+1); // # of cells in each section
						double unitcap= capRow1/param->numColSubArray;
						double unitres= resRow/param->numColSubArray;
						param->unitcap = unitcap;
						param->unitres = unitres;
						param->drivecapin = drivecapin;

						double capload_repeater1 = capRow1/(buffnum+1)+ drivecapin;
						double capload_repeater2 = capRow1/(buffnum+1);
						
						wlSwitchMatrix.CalculateLatency(1e20, capRow1/(buffnum+1)+ drivecapin, unitres * sectionnum, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);						
						
						double t1 =  targetdriveres * (drivecapout) * 0.69 
									+ unitcap * sectionnum * (0.69*targetdriveres + 0.38* unitres * sectionnum  )
									+ (unitcap * sectionnum * 0.69 + 0.69 *  targetdriveres  )* drivecapin;										
						double t2 =  targetdriveres * (drivecapout) * 0.69 
									+ unitcap * sectionnum * (0.69*targetdriveres + 0.38* unitres * sectionnum  );										
						double t3 = (drivecapout + drivecapin) * targetdriveres * 0.69;
						
						while (iterbuffnum >= 0){
							if (iterbuffnum  == 0) {
								bufferlatency += t2 + t3;
							}
							else {
								bufferlatency += t1+ t3;
							}
							iterbuffnum  = iterbuffnum  -1;
						}
					} else {
						bufferlatency =0;
						wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					}
					precharger.CalculateLatency(1e20, capCol, 1, numWriteOperationPerRow*numRow*activityRowWrite);
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}					
					if (numColMuxed > 1) {
						mux.CalculateLatency(1e20, 0, 1);						
						// 1.4 update: more arguments for muxdecoder.calculatelatency
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}					

					// Read

					// 1.4 update: SRAM column cap update (calibration with IMEC data)
					// 1.4 update: parallel mode discharge characteristics should be accounted
					// 1.4 update: SRAM row delay update

					// tau = (resCellAccess + resPullDown / m) * (capCellAccess + capCol) + resCol * capCol / 2;
					// m = enabled rows && weight 1

					// 1.4 update : new bitline model - needs check for Parallel mode
					double resPullDown = CalculateOnResistance(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech);
					double BLCap_perCell = capCol / numRow + capCellAccess; // Anni update
					double BLRes_perCell = resCol / numRow;
					double Elmore_BL = (resCellAccess + resPullDown) * BLCap_perCell * numRow   + BLCap_perCell * BLRes_perCell * numRow  * ( numRow +1 )  /2;

					colDelay = Elmore_BL * log(tech.vdd / (tech.vdd - cell.minSenseVoltage / 2));  

					if (CalculateclkFreq) {
						// 1.4 update - updated
						readLatency += MAX(wlSwitchMatrix.readLatency + bufferlatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) );
						readLatency += precharger.readLatency;
						// readLatency += colDelay;	   
						readLatency += multilevelSenseAmp.readLatency;

						            // --- INJECT THIS DEBUG BLOCK ---
                                        			//cout << "\n[DEBUG] SRAM Parallel Buffer Latency Breakdown:" << endl;
                                        			//cout << "  - wlDecoder: " << wlDecoder.readLatency << " s" << endl;
                                        			//cout << "  - precharger: " << precharger.readLatency << " s" << endl;
                                        			//cout << "  - colDelay: " << colDelay << " s" << endl;
                                        			//cout << "  - senseAmp: " << senseAmp.readLatency << " s" << endl;
                                        			//cout << "  - Total SRAM Parallel readLatency: " << readLatency << " s\n" << endl;
                                        			// -------------------------------


						param->rowdelay = wlSwitchMatrix.readLatency + bufferlatency;
						param->muxdelay = mux.readLatency+muxDecoder.readLatency;
						param->ADClatency = multilevelSenseAmp.readLatency;

						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: add partial sums
					if (numAdd > 1) {	
						adder.CalculateLatency(1e20, dff.capTgDrain, 1); // numRead = numColMuxed*(numAdd-1)
						dff.CalculateLatency(1e20, 1);	// numRead = numColMuxed*(numAdd-1)
					}
					if (numCellPerSynapse > 1) {
						shiftAddWeight.CalculateLatency(1); // numRead = (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse)
					}
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(1);	// numRead = ceil(numColMuxed/numCellPerSynapse)
					}					
					if (param->synchronous) {
						readLatencyADC = numColMuxed * numAdd;	// Anni update	
						// Anni update: hide readLatencyAccum by pipeline
						if (numAdd > 1) {	// adder is pipelined with ADC
							readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);	
						}	
						if (numCellPerSynapse > 1) {	// shiftAddWeight is pipelined with adder+ADC of one column
							readLatencyAccum += (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddWeight.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/numColMuxed, 0);	
						} 
						if (numReadPulse > 1) {	 // shiftAddInput is pipelined with adder+ADC+shiftaddweight of numCellPerSynapse columns
							readLatencyAccum += ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddInput.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/ceil(numColMuxed/numCellPerSynapse), 0);	
						} 
					} else {
						readLatencyADC = (precharger.readLatency + colDelay + multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency) * numColMuxed * (validated==true? param->beta : 1) * numAdd;
						readLatencyOther = MAX(wlSwitchMatrix.readLatency * numAdd, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) ) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed*(numAdd-1) + shiftAddWeight.adder.readLatency * (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) + \
											shiftAddInput.adder.readLatency * ceil(numColMuxed/numCellPerSynapse) - readLatencyADC - readLatencyOther, 0);
					}					
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
			} else if (BNNsequentialMode || XNORsequentialMode) {
				int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {

					// 1.4 update: new arguments for row decoder
					wlDecoder.CalculateLatency(1e20, capRow1, NULL, resRow, numCol, 1, numRow*activityRowWrite);
					precharger.CalculateLatency(1e20, capCol,1, numWriteOperationPerRow*numRow*activityRowWrite);				
					senseAmp.CalculateLatency(1);				
					
					// Read

					// 1.4 update : new bitline model - needs check for BNN/XNOR mode
					double resPullDown = CalculateOnResistance(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech);
					double BLCap_perCell = capCol / numRow + capCellAccess; // Anni update
					double BLRes_perCell = resCol / numRow;
					double Elmore_BL = (resCellAccess + resPullDown) * BLCap_perCell * numRow   + BLCap_perCell * BLRes_perCell * numRow  * ( numRow +1 )  /2;

					colDelay = Elmore_BL * log(tech.vdd / (tech.vdd - cell.minSenseVoltage / 2));  


					if (CalculateclkFreq) {
						readLatency += wlDecoder.readLatency;
						readLatency += precharger.readLatency;
						readLatency += colDelay;
						readLatency += senseAmp.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: hide readLatencyAccum by pipeline
					adder.CalculateLatency(1e20, dff.capTgDrain, 1);
					dff.CalculateLatency(1e20, 1);
					if (param->synchronous) {
						readLatencyADC = numReadOperationPerRow*numRow*activityRowRead;						
						// adder is pipelined with ADC
						readLatencyAccum = numReadOperationPerRow*(numRow*activityRowRead-1) * (ceil(adder.readLatency*clkFreq)-1);		
					} else {
						readLatencyADC = (precharger.readLatency + colDelay + senseAmp.readLatency) * numReadOperationPerRow*numRow*activityRowRead * (validated==true? param->beta : 1);						
						readLatencyOther = wlDecoder.readLatency * numRow*activityRowRead * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numReadOperationPerRow*(numRow*activityRowRead-1) - readLatencyADC - readLatencyOther, 0);
					}				
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}				
			} else if (BNNparallelMode || XNORparallelMode) {
				int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {
					wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					precharger.CalculateLatency(1e20, capCol, 1, numWriteOperationPerRow*numRow*activityRowWrite);					
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}
					if (numColMuxed > 1) {
						mux.CalculateLatency(1e20, 0, 1);

						// 1.4 update: new arguments for muxdecoder
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}				
					
					// Read
					// 1.4 update : new bitline model - needs check for BNN/XNOR parallel mode
					double resPullDown = CalculateOnResistance(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, inputParameter.temperature, tech);
					double BLCap_perCell = capCol / numRow + capCellAccess; // Anni update
					double BLRes_perCell = resCol / numRow;
					double Elmore_BL = (resCellAccess + resPullDown) * BLCap_perCell * numRow   + BLCap_perCell * BLRes_perCell * numRow  * ( numRow +1 )  /2;

					colDelay = Elmore_BL * log(tech.vdd / (tech.vdd - cell.minSenseVoltage / 2));
					
					if (CalculateclkFreq) {
						readLatency += MAX(wlSwitchMatrix.readLatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						readLatency += precharger.readLatency;
						// readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: add partial sums; hide readLatencyAccum by pipeline
					if (numAdd > 1) {	
						adder.CalculateLatency(1e20, dff.capTgDrain, 1);
						dff.CalculateLatency(1e20, 1);
					}
					if (param->synchronous) {
						readLatencyADC = numColMuxed * numAdd;
						// adder is pipelined with ADC
						readLatencyAccum = numColMuxed * (numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);
					} else {
						readLatencyADC = (precharger.readLatency + colDelay + multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency) * numColMuxed * (validated==true? param->beta : 1) * numAdd;
						readLatencyOther = MAX(wlSwitchMatrix.readLatency * numAdd, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0) ) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed * (numAdd-1) - readLatencyADC - readLatencyOther, 0);
					}
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}			
			}
	    } 
		else if (cell.memCellType == Type::Cap){ // nvCap added
			if (conventionalParallel) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double colRamp = 0;

				// 1.4 update: needs check - (capCol)*(cell.resMemCellAvg/(numRow));? 	Anni update: numRow->numRowParallel
				double tau = (capCol)*(cell.resMemCellAvg/(numRowParallel/2));

				// colDelay = horowitz(tau, 0, 1e20, &colRamp);
				// colDelay = tau * 0.2;  // assume the 15~20% voltage drop is enough for sensing
				// CAP update:
				colDelay = param->chargeDelay / 128 * numRowParallel;
				if (CalculateclkFreq || !param->synchronous) {				
					wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					if (numColMuxed>1) {
						mux.CalculateLatency(1e20, 0, 1);
						// 1.4 update
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}				
					if (CalculateclkFreq) {
						readLatency += MAX(wlSwitchMatrix.readLatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default		
						// cout<<"colDelay: "<<colDelay<<endl;			
						// cout<<"wlSwitchMatrix.readLatency: "<<wlSwitchMatrix.readLatency<<endl;	
						// cout<<"mux.readLatency: "<<mux.readLatency<<endl;		
						// cout<<"muxDecoder.readLatency: "<<muxDecoder.readLatency<<endl;	
						// cout<<"multilevelSenseAmp.readLatency: "<<multilevelSenseAmp.readLatency<<endl;	
						// cout<<"readLatency: "<<readLatency<<endl;		
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: add partial sums
					if (numAdd > 1) {	
						adder.CalculateLatency(1e20, dff.capTgDrain, 1); // numRead = numColMuxed*(numAdd-1)
						dff.CalculateLatency(1e20, 1);	// numRead = numColMuxed*(numAdd-1)
					}
					if (numCellPerSynapse > 1) {
						shiftAddWeight.CalculateLatency(1);	// numRead = (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse)
					}
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(1);	// numRead = ceil(numColMuxed/numCellPerSynapse)
					}
					if (param->synchronous) {
						readLatencyADC = numColMuxed * numAdd;	// Anni update	
						// Anni update: hide readLatencyAccum by pipeline
						if (numAdd > 1) {	// adder is pipelined with ADC
							readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);	
						}	
						if (numCellPerSynapse > 1) {	// shiftAddWeight is pipelined with adder+ADC of one column
							readLatencyAccum += (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddWeight.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/numColMuxed, 0);	
						} 
						if (numReadPulse > 1) {	 // shiftAddInput is pipelined with adder+ADC+shiftaddweight of numCellPerSynapse columns
							readLatencyAccum += ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddInput.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/ceil(numColMuxed/numCellPerSynapse), 0);	
						} 
					} else {
						readLatencyADC = (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency + colDelay) * numColMuxed * (validated==true? param->beta : 1) * numAdd;
						readLatencyOther = MAX((wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency) * numAdd, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0)) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed*(numAdd-1) + shiftAddWeight.adder.readLatency * (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) + \
											shiftAddInput.adder.readLatency * ceil(numColMuxed/numCellPerSynapse) - readLatencyADC - readLatencyOther, 0);
					}
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
			}		
		



	    } else if (cell.memCellType == Type::_2TnC) {

                        double bitsPerCell = param->bitsPerCell;
		        int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);

                        //Operation Counts (How many cycles to process the full array)
                        double numReadCells = (int)ceil((double)numCol/numColMuxed);   //In how many batches are you reading (8 cols share 1 ADC)
                        int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);  //hOW MANY WRITE CYCLES REQUIRED PER ROW

                        //Calculate Driver Latencies (RC Delay of the wires)
                        if (CalculateclkFreq || !param->synchronous) {

                                
				// Read Drivers

				// RSL (Col): Drives vertically down the rows
				//rslDecoder.CalculateLatency(1e20, capRSL, NULL, resCol, numRow, 1, 0);
				rslSwitchMatrix.CalculateLatency(1e20, capRSL, resCol, 1, 0);
				
				// RBL (Row): Drives horizontally across the columns
				//rblDecoder.CalculateLatency(1e20, capRBL, NULL, resRow, numCol, 1, 0);
				rblSwitchMatrix.CalculateLatency(1e20, capRBL, resRow, 1, 0);
				
				// Calculate Mux delay for Bitline Sensing
				if (numColMuxed > 1) {
				        mux.CalculateLatency(1e20, 0, 1);
				        muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
				}
				
				// WRITE PATH DRIVERS (WBL, SSL, WWL, WPL)
				
				// WBL (Row - Plane): Drives horizontally across the columns
				// wblDecoder.CalculateLatency(1e20, capWBL, NULL, resRow, numCol, 0, 1);
				// wblSwitchMatrix.CalculateLatency(1e20, capWBL, resRow, 1, 1);

				//wblDecoder.CalculateLatency(1e20, capWBL, NULL, resRow, numCol, 1, 1);
				wblPlaneDecoder.CalculateLatency(1e20, capWBL, NULL, resRow, numCol, 1, 1);
				// the two decode in parallel; the WBL path waits for the slower one
				wblDecoder.readLatency  = MAX(wblDecoder.readLatency,  wblPlaneDecoder.readLatency);
				wblDecoder.writeLatency = MAX(wblDecoder.writeLatency, wblPlaneDecoder.writeLatency);
				//wblSwitchMatrix.CalculateLatency(1e20, capWBL, resRow, 1, 1);
				//wblDecoder.CalculateLatency(1e20, capWBL, NULL, resWBL, numCol, 1, 1);
				wblSwitchMatrix.CalculateLatency(1e20, capWBL, resWBL, 1, 1);

				// SSL (Col): Drives vertically down the rows
				//sslDecoder.CalculateLatency(1e20, capSSL, NULL, resCol, numRow, 0, 1);
				sslSwitchMatrix.CalculateLatency(1e20, capSSL, resCol, 1, 1);
				
				// WWL (Row): Drives horizontally across the columns
				//wwlDecoder.CalculateLatency(1e20, capWWL, NULL, resRow, numCol, 0, 1);
				wwlSwitchMatrix.CalculateLatency(1e20, capWWL, resRow, 1, 1);
				
				// WPL (Row): Drives horizontally across the columns
				//wplDecoder.CalculateLatency(1e20, capWPL, NULL, resRow, numCol, 0, 1);
				wplSwitchMatrix.CalculateLatency(1e20, capWPL, resRow, 1, 1);
				
				
				double capDecOut = 4.0 * gateCap;      // ~0.1 fF, one TG select input 

				rslDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numRow, 1, 0);
				rblDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numCol, 1, 0);
				wblDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numCol, 1, 1);
				sslDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numRow, 0, 1);
				wwlDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numCol, 0, 1);
				wplDecoder.CalculateLatency(1e20, capDecOut, NULL, 0, numCol, 0, 1);
				
				//Sensisng Circuits
                                if (param->SARADC) {
                                        sarADC.CalculateLatency(1);
                                } else {
					multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
                                        multilevelSAEncoder.CalculateLatency(1e20, 1);

					multilevelSenseAmp.readLatency = 1e-9;
                                }

				double Trcd = 85e-9;
				double Trp = 80e-9;

                                //Total Read Latency
                                if (CalculateclkFreq) {

					colDelay = param->chargeDelay;
        				readLatency += colDelay;

					// QNDRO SELF-RESTORE PENALTY 
					if (cell.mem_rdo == Type::qndro) {

	                                    if (writeBack) {

						    // 1. Calculate Single-Row Read-Out Latency
	                                	    // Parameters: (columnResistance, numColMuxed, numRead)
	                                	    // numColMuxed = 1: Dedicated SA per column (no multiplexing)
	                                	    // numRead = 1: Single sensing operation
	                                	    currentSenseAmp.CalculateLatency(columnResistance, 1, 1);
	
	                                	    // 2. Scale sensing latency for the entire subarray (row-by-row)
	                                	    double senseLatencyTotal = currentSenseAmp.readLatency * numRow * (int)param->bitsPerCell;
	
	                                	        // Find the max latency of the write drivers
	                                	    double maxWriteDriverLatency = 0;
	                                	    maxWriteDriverLatency = MAX(maxWriteDriverLatency, wblSwitchMatrix.writeLatency + wblDecoder.writeLatency);
	                                	    maxWriteDriverLatency = MAX(maxWriteDriverLatency, sslSwitchMatrix.writeLatency + sslDecoder.writeLatency);
	                                	    maxWriteDriverLatency = MAX(maxWriteDriverLatency, wwlSwitchMatrix.writeLatency + wwlDecoder.writeLatency);
	                                	    maxWriteDriverLatency = MAX(maxWriteDriverLatency, wplSwitchMatrix.writeLatency + wplDecoder.writeLatency);
	
	                                	    double writePulseTime = 10e-9;
							//double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
							//double writePulseTime = param->nlsTauInf
							//        * exp(pow(param->nlsEa / eWrField, param->nlsAlpha));
						    int numWriteOperationPerRow = (int)ceil((double)numCol/numWriteCellPerOperationNeuro);
	
	                                	    // A full array write-back requires rewriting every row in the subarray sequentially
	                                	    // double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * bitsPerCell * numRow;

						    double numRowRestored = numRow * activityRowRead;
						    double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * numRowRestored;
						    restoreLatencyCore       =  writePulseTime                          * numRowRestored;
	
						    // Total Refresh Latency = Read-Out (Sensing) + Write Phase
	                                	    double writeBackLatency = senseLatencyTotal + writePhaseLatency;
	
	                                	    readLatency += writeBackLatency;
	                                	}

					}

					else if (cell.mem_rdo == Type::dro) {
    					    // 1. Calculate Single-Row Read-Out Latency
                                            // Parameters: (columnResistance, numColMuxed, numRead)
                                            // numColMuxed = 1: Dedicated SA per column (no multiplexing)
                                            // numRead = 1: Single sensing operation
                                            currentSenseAmp.CalculateLatency(columnResistance, 1, 1);

                                            // 2. Scale sensing latency for the entire subarray (row-by-row)
                                            double senseLatencyTotal = currentSenseAmp.readLatency * numRow * bitsPerCell;

						
					    // Find the max latency of the write drivers needed for write-back
    					    double maxWriteDriverLatency = 0;
    					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, wblSwitchMatrix.writeLatency + wblDecoder.writeLatency);
    					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, sslSwitchMatrix.writeLatency + sslDecoder.writeLatency);
    					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, wwlSwitchMatrix.writeLatency + wwlDecoder.writeLatency);
    					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, wplSwitchMatrix.writeLatency + wplDecoder.writeLatency);

    					    maxWriteDriverLatency = MAX(maxWriteDriverLatency, wblSwitchMatrix.readLatency + wblDecoder.readLatency);
                                            maxWriteDriverLatency = MAX(maxWriteDriverLatency, sslSwitchMatrix.readLatency + sslDecoder.readLatency);
                                            maxWriteDriverLatency = MAX(maxWriteDriverLatency, wwlSwitchMatrix.readLatency + wwlDecoder.readLatency);
                                            maxWriteDriverLatency = MAX(maxWriteDriverLatency, wplSwitchMatrix.readLatency + wplDecoder.readLatency);

					    
					    double writePulseTime = 10e-9;
				        //     double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
				        //     double writePulseTime = param->nlsTauInf
				        // 			* exp(pow(param->nlsEa / eWrField, param->nlsAlpha));
					    int numWriteOperationPerRow = (int)ceil((double)numCol/numWriteCellPerOperationNeuro);

                                            // A full array write-back requires rewriting every row in the subarray sequentially
                                            //double writePhaseLatency = (maxWriteDriverLatency + writePulseTime + Twr) * numRow ;
                                            // double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * numRow * bitsPerCell;
					    
				            double numRowRestored = numRow * activityRowRead;
                                            double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * numRowRestored;
                                            restoreLatencyCore       =  writePulseTime                          * numRowRestored;

					    // Total Refresh Latency = Read-Out (Sensing) + Write Phase
                                            double writeBackLatency = senseLatencyTotal + writePhaseLatency;

                                            readLatency += writeBackLatency;
    					}


					// Map 3D line parasitics to NeuroSim's calculated array parasitics
					double resReadTransistor = cell.resMemCellOn; // Assuming worst-case ON resistance

                                        // 1. Setup Delay: WBL and RBL are driven simultaneously
        				// Calculate RC delay for WBL (Input line)
        				double delay_WBL = 0.69 * resRow * capWBL; 

        				// Calculate RC delay for RBL (Bias line)
        				double delay_RBL = 0.69 * resRow * capRBL;

        				// The setup time is governed by the slower of the two lines
        				double delay_setup = max(delay_WBL, delay_RBL);

        				// 2. Sensing Delay: Time taken for RSL to develop detectable current/voltage
        				// This depends on the TR on-resistance and RSL capacitance
        				//double delay_RSL = 0.69 * (resReadTransistor + resCellRSL * numRow) * capRSL;

        				// 3. Peripheral Delay
        				//double delay_peripherals = wlDecoder.readLatency + mux.readLatency + multilevelSenseAmp.readLatency;

        				// Total Read Latency
        				//readLatency = delay_setup + delay_RSL + delay_peripherals;
        				readLatency += delay_setup;
					lineSetupCore    = delay_setup;
					// readLatencyCore  = colDelay + delay_setup + restoreLatencyCore;
					{
					    double restoreDutyL = 0.0;
					    if      (cell.mem_rdo == Type::dro)   restoreDutyL = 1.0;
					    else if (cell.mem_rdo == Type::qndro) restoreDutyL = 1.0 / param->qndroRefreshInterval;
					    restoreLatencyCore = restoreDutyL * 10e-9 * numRow * activityRowRead;
					}
					readLatencyCore  = colDelay + delay_setup + restoreLatencyCore;


						currentSenseAmp.CalculateLatency(columnResistance, 1, 1);
						senseLatencyCore = currentSenseAmp.readLatency;
						readLatencySense = readLatencyCore + senseLatencyCore;

					
					//Driver Delay: Time to activate RSL
                                        // 2TnC Driver Delay: WBL (Input), RBL (Drain Bias), RSL (Source Bias)
					double wblPath = wblDecoder.readLatency + wblSwitchMatrix.readLatency;
					double rblPath = rblDecoder.readLatency + rblSwitchMatrix.readLatency;
					double rslPath = rslDecoder.readLatency + rslSwitchMatrix.readLatency;
					
					double driverDelay = MAX(wblPath, MAX(rblPath, rslPath))
					                     + ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);
					
					readLatency += driverDelay;

                                        //Sensing Delay (ADC + Encoder)
                                        if (param->SARADC) {
                                                readLatency += sarADC.readLatency;
                                        } else {
                                                readLatency += multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency;
                                        }

                                        //Scaling factor (beta)
                                        readLatency *= (validated==true? param->beta : 1);

                                }
			}


                        //Write Latency
                        if (!CalculateclkFreq) {
                                //Max of all 4 switch matrices
                                double maxDriverLatency = 0;
                                maxDriverLatency = MAX(maxDriverLatency, wblSwitchMatrix.writeLatency + wblDecoder.writeLatency);
                                maxDriverLatency = MAX(maxDriverLatency, sslSwitchMatrix.writeLatency + sslDecoder.writeLatency);
                                maxDriverLatency = MAX(maxDriverLatency, wwlSwitchMatrix.writeLatency + wwlDecoder.writeLatency);
                                maxDriverLatency = MAX(maxDriverLatency, wplSwitchMatrix.writeLatency + wplDecoder.writeLatency);

                                //Calculate Capacitor Charging Time (RC Delay)
                                //tau = R * C
                                //double totalCap = cell.numCapacitors * cell.capacitance; //Load of the NvCAPs
                                //double tau = cell.resistanceAccess * totalCap;           //Time constant

                                //double writePulseTime = 3.0 * tau;

				//double writePulseTime = cell.writePulseWidth;
				double writePulseTime = 10e-9;

				//double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
				//double writePulseTime = param->nlsTauInf
				//        * exp(pow(param->nlsEa / eWrField, param->nlsAlpha));


				// 2-Cycle Parallel Write Scheme
				// Cycle 1: Write all '1's (WPL = High, Target WBLs = Low, Inhibit WBLs = High)
				// Cycle 2: Write all '0's (WPL = Low, Target WBLs = High, Inhibit WBLs = Low)
				// int writeCyclesPerCell = 2;
				
				// Total Write Latency = (Driver Setup + Pulse Duration) * 2 cycles * Num Column Operations
				//writeLatency = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow * numRow;

				double numRowWritten = numRow * activityRowWrite;
				writeLatency     = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow * numRowWritten;
				writeLatencyCore =  writePulseTime                     * numWriteOperationPerRow * numRowWritten;

                                // writeLatency = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow;


                                //Read Latency (Non-Clocked Mode)
    				if (param->synchronous) {
    				        readLatencyADC = numColMuxed;
    				        readLatencyAccum = 0;
    				        if (numAdd > 1) {
    				                readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);
    				        }
    				} else {
    				        double sensingDelay = (param->SARADC ? sarADC.readLatency : (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency));

					// 2TnC Driver Delay: WBL (Input), RBL (Drain Bias), RSL (Source Bias)
					double wblPath = wblDecoder.readLatency + wblSwitchMatrix.readLatency;
					double rblPath = rblDecoder.readLatency + rblSwitchMatrix.readLatency;
					double rslPath = rslDecoder.readLatency + rslSwitchMatrix.readLatency;
					
					double driverDelay = MAX(wblPath, MAX(rblPath, rslPath))
					                     + ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);

    				        readLatencyADC = (driverDelay + sensingDelay) * numColMuxed * (validated==true? param->beta : 1);
    				        readLatencyAccum = MAX(((numAdd>1) ? adder.readLatency * numColMuxed*(numAdd-1) : 0) 
						+ ((numCellPerSynapse > 1) ? shiftAddWeight.readLatency : 0)
						+ ((numReadPulse > 1) ? shiftAddInput.readLatency : 0) - readLatencyADC, 0.0);
    				}

    				//Accumulation Latency (Digital Adder/Shift)
    				if (numAdd > 1) {
    				        adder.CalculateLatency(1e20, dff.capTgDrain, 1);
    				        dff.CalculateLatency(1e20, 1);
    				}
    				if (numCellPerSynapse > 1) 
					shiftAddWeight.CalculateLatency(numColMuxed);
    				if (numReadPulse > 1) 
					shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));

    				readLatencyOther = param->chargeDelay + lineSetupCore;
    				readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				// readLatency += restoreLatencyCore;   /* REVERTED: ProcessingUnit adds the
				//                                        write-back itself via WRITEBACKCYCLE.
				//                                        Touching readLatency here double-counts
				//                                        it for CIM. */
			}

			 // cout << "2TnC Read Latency Breakdown " << endl;
                         // cout << "numReadCellPerOperationNeuro: " << numReadCellPerOperationNeuro << endl;
                         // cout << "numReadOperationPerRow: " << numReadOperationPerRow << endl;
                         // cout << "numReadCells: " << numReadCells << endl;
                         // cout << "rslDecoder.readLatency: " << rslDecoder.readLatency << endl;
                         // cout << "rslSwitchMatrix.readLatency: " << rslSwitchMatrix.readLatency << endl;
                         // cout << "mux.readLatency: " << mux.readLatency << endl;
                         // cout << "muxDecoder.readLatency: " << muxDecoder.readLatency << endl;
                         // cout << "sarADC.readLatency: " << sarADC.readLatency << endl;
                         // // cout << "multilevelSenseAmp.readLatency: " << multilevelSenseAmp.readLatency << endl;
                         // // cout << "multilevelSAEncoder.readLatency: " << multilevelSAEncoder.readLatency << endl;
                         // cout << "readLatencyADC: " << readLatencyADC << endl;
                         // cout << "readLatencyAccum: " << readLatencyAccum << endl;
                         // 
			 // cout << "2TnC Write Latency Breakdown " << endl;
			 // cout << "CapWBL: " << capWBL << endl;
			 // cout << "CapWWL: " << capWWL << endl;
			 // cout << "CapWPL: " << capWPL << endl;
			 // cout << "CapSSL: " << capSSL << endl;
			 // cout << "wblSwitchMatrix.writeLatency: " << wblSwitchMatrix.writeLatency << endl;
			 // cout << "wblSwitchMatrix.writeLatency: " << wblSwitchMatrix.writeLatency << endl;
			 // cout << "wblSwitchMatrix.writeLatency: " << wblSwitchMatrix.writeLatency << endl;
			 // cout << "wplSwitchMatrix.writeLatency: " << wplSwitchMatrix.writeLatency << endl;
			 // cout << "wwlSwitchMatrix.writeLatency: " << wwlSwitchMatrix.writeLatency << endl;
			 // cout << "sslSwitchMatrix.writeLatency: " << sslSwitchMatrix.writeLatency << endl;
			 // cout << "wblDecoder.writeLatency: " << wblDecoder.writeLatency << endl;
			 // cout << "wwlDecoder.writeLatency: " << wwlDecoder.writeLatency << endl;
			 // cout << "wplDecoder.writeLatency: " << wplDecoder.writeLatency << endl;
			 // cout << "sslDecoder.writeLatency: " << sslDecoder.writeLatency << endl;
			 // cout << "numWriteOperationPerRow: " << numWriteOperationPerRow << endl;
			 // //cout << ": " <<  << endl;
		
		
		} else if (cell.memCellType == Type::_1TnC) {
        		double bitsPerCell = param->bitsPerCell;
			int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);

                        //Operation Counts (How many cycles to process the full array)
                        double numReadCells = (int)ceil((double)numCol/numColMuxed);   //In how many batches are you reading (8 cols share 1 ADC)
                        int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);  //hOW MANY WRITE CYCLES REQUIRED PER ROW

                        //Calculate Driver Latencies (RC Delay of the wires)
                        if (CalculateclkFreq || !param->synchronous) {

				// READ DRIVERS (BL)
				// BL (Row - Plane): Drives horizontally across the columns
				blDecoder.CalculateLatency(1e20, capBL, NULL, resRow, numCol, 0, 1);
				blSwitchMatrix.CalculateLatency(1e20, capBL, resRow, 0, 1);
				
				// Calculate Mux delay
				if (numColMuxed > 1) {
				        mux.CalculateLatency(1e20, 0, 1);
				        muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
				}
				
				// WRITE PATH DRIVERS (WL, PL)
				// 2 lines are driven. The slowest one gates the operation.
				// WL (Row): Wordline (Transistor Gate) - Drives horizontally across the columns
				wlDecoder.CalculateLatency(1e20, capWL, NULL, resRow, numCol, 1, 1);
				wlSwitchMatrix.CalculateLatency(1e20, capWL, resRow, 1, 1);
				
				// PL (Column): Write Plate Line - Drives vertically down the rows
				plDecoder.CalculateLatency(1e20, capPL, NULL, resCol, numRow, 1, 1);
				plSwitchMatrix.CalculateLatency(1e20, capPL, resCol, 1, 1);
				
				
                                //Sensisng Circuits
                                if (param->SARADC) {
                                        sarADC.CalculateLatency(1);
                                } else {
                                        multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
                                        multilevelSAEncoder.CalculateLatency(1e20, 1);

                                        multilevelSenseAmp.readLatency = 1e-9;
                                }

                                double Trc = 185e-9;
                                double Trcd = 85e-9;
                                double Trp = 80e-9;

                                //Total Read Latency
                                if (CalculateclkFreq) {

                                        colDelay = param->chargeDelay;
                                        readLatency += colDelay;

                                        // QNDRO SELF-RESTORE PENALTY 
                                        if (cell.mem_rdo == Type::qndro) {

                                            if (writeBack) {

                                                    // 1. Calculate Single-Row Read-Out Latency
                                                    // Parameters: (columnResistance, numColMuxed, numRead)
                                                    // numColMuxed = 1: Dedicated SA per column (no multiplexing)
                                                    // numRead = 1: Single sensing operation
                                                    currentSenseAmp.CalculateLatency(columnResistance, 1, 1);

                                                    // 2. Scale sensing latency for the entire subarray (row-by-row)
                                                    double senseLatencyTotal = currentSenseAmp.readLatency * numRow * (int)param->bitsPerCell;

                                                        // Find the max latency of the write drivers
                                                    double maxWriteDriverLatency = 0;
                                                    maxWriteDriverLatency = MAX(maxWriteDriverLatency, blSwitchMatrix.writeLatency + blDecoder.writeLatency);
                                                    maxWriteDriverLatency = MAX(maxWriteDriverLatency, wlSwitchMatrix.writeLatency + wlDecoder.writeLatency);
                                                    maxWriteDriverLatency = MAX(maxWriteDriverLatency, plSwitchMatrix.writeLatency + plDecoder.writeLatency);

                                                    double writePulseTime = 10e-9;
							// double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
							// double writePulseTime = param->nlsTauInf
							//         * exp(pow(param->nlsEa / eWrField, param->nlsAlpha));

						    int numWriteOperationPerRow = (int)ceil((double)numCol/numWriteCellPerOperationNeuro);

                                                    // A full array write-back requires rewriting every row in the subarray sequentially
                                                    // double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * bitsPerCell * numRow;

						    	double numRowRestored = numRow * activityRowRead;
							double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * numRowRestored;
							restoreLatencyCore       =  writePulseTime                          * numRowRestored;

                                                    // Total Refresh Latency = Read-Out (Sensing) + Write Phase
                                                    double writeBackLatency = senseLatencyTotal + writePhaseLatency;

                                                    readLatency += writeBackLatency;
                                                }

                                        }

                                        else if (cell.mem_rdo == Type::dro) {
                                            // // 1. Calculate Single-Row Read-Out Latency
                                            // // Parameters: (columnResistance, numColMuxed, numRead)
                                            // // numColMuxed = 1: Dedicated SA per column (no multiplexing)
                                            // // numRead = 1: Single sensing operation
                                            // currentSenseAmp.CalculateLatency(columnResistance, 1, 1);

                                            // // 2. Scale sensing latency for the entire subarray (row-by-row)
                                            // double senseLatencyTotal = currentSenseAmp.readLatency * numRow * bitsPerCell;


                                            // // Find the max latency of the write drivers needed for write-back
                                            // double maxWriteDriverLatency = 0;
                                            // // maxWriteDriverLatency = MAX(maxWriteDriverLatency, blSwitchMatrix.writeLatency + blDecoder.writeLatency);
                                            // // maxWriteDriverLatency = MAX(maxWriteDriverLatency, wlSwitchMatrix.writeLatency + wlDecoder.writeLatency);
                                            // // maxWriteDriverLatency = MAX(maxWriteDriverLatency, plSwitchMatrix.writeLatency + plDecoder.writeLatency);

                                            // maxWriteDriverLatency = MAX(maxWriteDriverLatency, blSwitchMatrix.readLatency + blDecoder.readLatency);
					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, wlSwitchMatrix.readLatency + wlDecoder.readLatency);
					    // maxWriteDriverLatency = MAX(maxWriteDriverLatency, plSwitchMatrix.readLatency + plDecoder.readLatency);

					    // double writePulseTime = 10e-9;
                                            // int numWriteOperationPerRow = (int)ceil((double)numCol/numWriteCellPerOperationNeuro);

                                            // // A full array write-back requires rewriting every row in the subarray sequentially
                                            // //double writePhaseLatency = (maxWriteDriverLatency + writePulseTime + Twr) * numRow ;
                                            // double writePhaseLatency = (maxWriteDriverLatency + writePulseTime) * numRow * bitsPerCell;

                                            // // Total Refresh Latency = Read-Out (Sensing) + Write Phase
                                            // double writeBackLatency = senseLatencyTotal + writePhaseLatency;

                                            //double writeBackLatency = Trc * (numRow * activityRowRead);
                                             double writeBackLatency = Trc ;

					  //  cout << "activityRowRead: " << activityRowRead << endl;
					  //  cout << "writeBackLatency: " << writeBackLatency << endl;

					    readLatency += writeBackLatency;
                                        }


					{
					    double restoreDutyL = 0.0;
					    if      (cell.mem_rdo == Type::dro)   restoreDutyL = 1.0;
					    else if (cell.mem_rdo == Type::qndro) restoreDutyL = 1.0 / param->qndroRefreshInterval;
					    restoreLatencyCore = restoreDutyL * 10e-9 * numRow * activityRowRead;
					}

					double delay_BL    = 0.69 * resRow * capBL;
					double delay_WL    = 0.69 * resRow * capWL;
					double delay_setup = max(delay_BL, delay_WL);

					readLatency     += delay_setup;
					lineSetupCore    = delay_setup;
					readLatencyCore  = colDelay + delay_setup + restoreLatencyCore;

					currentSenseAmp.CalculateLatency(columnResistance, 1, 1);
					senseLatencyCore = currentSenseAmp.readLatency;
					readLatencySense = readLatencyCore + senseLatencyCore;

					double blPathDelay = blDecoder.readLatency + blSwitchMatrix.readLatency;
					double wlPathDelay = wlDecoder.readLatency + wlSwitchMatrix.readLatency;
					double plPathDelay = plDecoder.readLatency + plSwitchMatrix.readLatency; 
					
					// Add the slowest parallel driver path to the total latency
					double driverDelay = MAX(blPathDelay, wlPathDelay);
					driverDelay = MAX(driverDelay, plPathDelay);

					readLatency += driverDelay;


                                        //Mux Delay
                                        readLatency += ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);

                                        //Sensing Delay (ADC + Encoder)
                                        if (param->SARADC) {
                                                readLatency += sarADC.readLatency;
                                        } else {
                                                readLatency += multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency;
                                        }

                                        //Scaling factor (beta)
                                        readLatency *= (validated==true? param->beta : 1);

                                }
			}


                        //Write Latency
                        if (!CalculateclkFreq) {
                                double Trc = 185e-9;
				{

				
				//Max of all 3 switch matrices
                                double maxDriverLatency = 0;
                                maxDriverLatency = MAX(maxDriverLatency, blSwitchMatrix.writeLatency + blDecoder.writeLatency);
                                maxDriverLatency = MAX(maxDriverLatency, wlSwitchMatrix.writeLatency + wlDecoder.writeLatency);
                                maxDriverLatency = MAX(maxDriverLatency, plSwitchMatrix.writeLatency + plDecoder.writeLatency);

                                //Calculate Capacitor Charging Time (RC Delay)
                                //tau = R * C
                                //double totalCap = cell.numCapacitors * cell.capacitance; //Load of the NvCAPs
                                //double tau = cell.resistanceAccess * totalCap;           //Time constant

                                //double writePulseTime = 3.0 * tau;

                                //double writePulseTime = cell.writePulseWidth;
                                double writePulseTime = 10e-9;
			//	double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
			//	double writePulseTime = param->nlsTauInf
			//	        * exp(pow(param->nlsEa / eWrField, param->nlsAlpha));

                                // 2-Cycle Parallel Write Scheme
                                // Cycle 1: Write all '1's (WPL = High, Target WBLs = Low, Inhibit WBLs = High)
                                // Cycle 2: Write all '0's (WPL = Low, Target WBLs = High, Inhibit WBLs = Low)
                                // int writeCyclesPerCell = 2;

                                // Total Write Latency = (Driver Setup + Pulse Duration) * 2 cycles * Num Column Operations
                                //writeLatency = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow * numRow;

				double numRowWritten = numRow * activityRowWrite;
                                writeLatency     = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow * numRowWritten;
                                writeLatencyCore =  writePulseTime                     * numWriteOperationPerRow * numRowWritten;


                                // writeLatency = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow;




				//Read Latency (Non-Clocked Mode)
                                if (param->synchronous) {
                                        readLatencyADC = numColMuxed;
                                        readLatencyAccum = 0;
                                        if (numAdd > 1) {
                                                readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);
                                        }
                                } else {
                                        double sensingDelay = (param->SARADC ? sarADC.readLatency : (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency));

                                        double blPathDelay = blDecoder.readLatency + blSwitchMatrix.readLatency;
                                        double wlPathDelay = wlDecoder.readLatency + wlSwitchMatrix.readLatency;
                                        double plPathDelay = plDecoder.readLatency + plSwitchMatrix.readLatency;

                                        // Add the slowest parallel driver path to the total latency
                                        double driverDelay = MAX(blPathDelay, wlPathDelay);
                                        driverDelay = MAX(driverDelay, plPathDelay);

					// Add the slowest parallel driver path to the total latency
                                        driverDelay += ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);

                                        readLatencyADC = (driverDelay + sensingDelay) * numColMuxed * (validated==true? param->beta : 1);
                                        readLatencyAccum = MAX(((numAdd>1) ? adder.readLatency * numColMuxed*(numAdd-1) : 0)
                                                + ((numCellPerSynapse > 1) ? shiftAddWeight.readLatency : 0)
                                                + ((numReadPulse > 1) ? shiftAddInput.readLatency : 0) - readLatencyADC, 0.0);
                                }

                                //Accumulation Latency (Digital Adder/Shift)
                                if (numAdd > 1) {
                                        adder.CalculateLatency(1e20, dff.capTgDrain, 1);
                                        dff.CalculateLatency(1e20, 1);
                                }
                                if (numCellPerSynapse > 1)
                                        shiftAddWeight.CalculateLatency(numColMuxed);
                                if (numReadPulse > 1)
                                        shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));

                                readLatencyOther = param->chargeDelay + lineSetupCore;
                                readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
                                if (cell.mem_rdo == Type::dro) readLatency += Trc;   /* destructive read-out penalty */
				}
                        }

			// cout << "\n================ 1TnC Block Latency ================" << endl;
        		// cout << "--- Read Path ---" << endl;
        		// if (wlDecoder.initialized) cout << "WL Decoder Read Latency:       " << wlDecoder.readLatency << " ns" << endl;
        		// if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Read Latency: " << wlSwitchMatrix.readLatency << " ns" << endl;
        		// if (blDecoder.initialized) cout << "BL Decoder Read Latency:       " << blDecoder.readLatency << " ns" << endl;
        		// if (blSwitchMatrix.initialized) cout << "BL Switch Matrix Read Latency: " << blSwitchMatrix.readLatency << " ns" << endl;
        		// if (plDecoder.initialized) cout << "PL Decoder Read Latency:       " << plDecoder.readLatency << " ns" << endl;
        		// if (plSwitchMatrix.initialized) cout << "PL Switch Matrix Read Latency: " << plSwitchMatrix.readLatency << " ns" << endl;

        		// if (numColMuxed > 1) {
        		//     if (muxDecoder.initialized) cout << "MUX Decoder Read Latency:      " << muxDecoder.readLatency << " ns" << endl;
        		//     if (mux.initialized) cout << "MUX Read Latency:              " << mux.readLatency << " ns" << endl;
        		// }
        		// if (currentSenseAmp.initialized) {
        		//     cout << "Current Sense Amp Read Latency:" << currentSenseAmp.readLatency << " ns" << endl;
        		// }
        		// if (param->SARADC) {
        		//     if (sarADC.initialized) cout << "SAR ADC Read Latency:          " << sarADC.readLatency << " ns" << endl;
        		// } else {
        		//     if (multilevelSenseAmp.initialized) cout << "Multilevel SA Read Latency:    " << multilevelSenseAmp.readLatency << " ns" << endl;
        		//     if (multilevelSAEncoder.initialized) cout << "Multilevel SA Encoder Latency: " << multilevelSAEncoder.readLatency << " ns" << endl;
        		// }
        		// if (numAdd > 1) {
        		//     if (adder.initialized) cout << "Adder Read Latency:            " << adder.readLatency << " ns" << endl;
        		//     if (dff.initialized) cout << "DFF Read Latency:              " << dff.readLatency << " ns" << endl;
        		// }
        		// if (numCellPerSynapse > 1) {
        		//     if (shiftAddWeight.initialized) cout << "ShiftAdd Weight Read Latency:  " << shiftAddWeight.readLatency << " ns" << endl;
        		// }
        		// if (numReadPulse > 1) {
        		//     if (shiftAddInput.initialized) cout << "ShiftAdd Input Read Latency:   " << shiftAddInput.readLatency << " ns" << endl;
        		// }

        		// cout << "--- Write Path ---" << endl;
        		// if (wlDecoder.initialized) cout << "WL Decoder Write Latency:      " << wlDecoder.writeLatency << " ns" << endl;
        		// if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Write Latency:" << wlSwitchMatrix.writeLatency << " ns" << endl;
        		// if (blDecoder.initialized) cout << "BL Decoder Write Latency:      " << blDecoder.writeLatency << " ns" << endl;
        		// if (blSwitchMatrix.initialized) cout << "BL Switch Matrix Write Latency:" << blSwitchMatrix.writeLatency << " ns" << endl;
        		// if (plDecoder.initialized) cout << "PL Decoder Write Latency:      " << plDecoder.writeLatency << " ns" << endl;
        		// if (plSwitchMatrix.initialized) cout << "PL Switch Matrix Write Latency:" << plSwitchMatrix.writeLatency << " ns" << endl;

        		// if (cell.writeVoltage > 1.5) {
        		//     if (wlLevelShifter.initialized) cout << "WL Level Shifter Write Latency:" << wlLevelShifter.writeLatency << " ns" << endl;
        		//     if (blLevelShifter.initialized) cout << "BL Level Shifter Write Latency:" << blLevelShifter.writeLatency << " ns" << endl;
        		//     if (plLevelShifter.initialized) cout << "PL Level Shifter Write Latency:" << plLevelShifter.writeLatency << " ns" << endl;
        		// }

        		// // --- Destructive Read-Out Write-Back Check ---
        		// if (cell.mem_rdo == Type::dro) {
        		//     cout << "---------------- DRO Penalty -----------------------" << endl;
        		//     if (wlDecoder.initialized) cout << "[Write-Back] WL Decoder:       " << wlDecoder.writeLatency << " ns" << endl;
        		//     if (wlSwitchMatrix.initialized) cout << "[Write-Back] WL SwitchMatrix:  " << wlSwitchMatrix.writeLatency << " ns" << endl;
        		//     if (blDecoder.initialized) cout << "[Write-Back] BL Decoder:       " << blDecoder.writeLatency << " ns" << endl;
        		//     if (blSwitchMatrix.initialized) cout << "[Write-Back] BL SwitchMatrix:  " << blSwitchMatrix.writeLatency << " ns" << endl;
        		//     if (plDecoder.initialized) cout << "[Write-Back] PL Decoder:       " << plDecoder.writeLatency << " ns" << endl;
        		//     if (plSwitchMatrix.initialized) cout << "[Write-Back] PL SwitchMatrix:  " << plSwitchMatrix.writeLatency << " ns" << endl;
        		//     cout << "----------------------------------------------------" << endl;
        		// } else {
        		//     cout << "Mode: NDRO (No Write-Back Latency Penalty)" << endl;
        		// }
        		// cout << "====================================================" << endl;


		} else if (cell.memCellType == Type::_1T1C) {
        		 int numReadOperationPerRow = (int)ceil((double)numCol/numReadCellPerOperationNeuro);
        		 int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);

        		 double Trc = 60e-9;

			 if (CalculateclkFreq || !param->synchronous) {
        		     // BL (Col Mode)
        		     blDecoder.CalculateLatency(1e20, capBL, NULL, resCol, numRow, 1, 0);
        		     blSwitchMatrix.CalculateLatency(1e20, capBL, resCol, 1, 0);

        		     if (numColMuxed > 1) {
        		         mux.CalculateLatency(1e20, 0, 1);
        		         muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
        		     }

        		     // WL (Row Mode)
        		     wlDecoder.CalculateLatency(1e20, capWL, NULL, resRow, numCol, 0, 1);
        		     wlSwitchMatrix.CalculateLatency(1e20, capWL, resRow, 1, 1);

        		     if (param->SARADC) {
        		         sarADC.CalculateLatency(1);
        		     } else {
        		         multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
        		         multilevelSAEncoder.CalculateLatency(1e20, 1);
        		     }

        		     if (CalculateclkFreq) {

				     } else if (cell.memCellType == Type::_1T1C) {
						// 1. Normal Inference Operation (Destructive Read-Out)

						// Row Cycle Time (Trc): Total time to access ONE row and write it back
						double Trc = 60e-9;

						// Split the cycle into the Read phase and the Write-Back phase
						double arrayReadPhase = Trc / 2.0;
						double writeBackPhase = Trc ;

						// CIM parallel read allows simultaneous charge dumping,
						// but the Write-Back MUST be completely serialized row-by-row.
						double serializedRows = numRow * activityRowRead;

						// Total DRO penalty is the write-back time * every row that was destroyed
						double totalWriteBackLatency = writeBackPhase * serializedRows;

						// Add the explicit serialized write-back hardware lockout time
						readLatency += totalWriteBackLatency;

						// 2. Distributed Refresh Latency Penalty (Statistical Availability)

						// Retention time before 1T1C loses charge
						// double tRefresh = 32e-3; // 32 ms
						double tRefresh = 20.48e-3; // 20.48 ms

						// Calculate how often a single row must be refreshed
						double tRefreshPerRow = tRefresh / numRow;

						// A refresh event only restores ONE row, so it only takes 1 Trc
						double refreshEventLatency = Trc;

						// Calculate Refresh Duty Cycle (Hardware Lockout Percentage)
						double refreshDutyCycle = refreshEventLatency / tRefreshPerRow;

						// Apply the Statistical Penalty to the Total Read Latency
						readLatency = readLatency * (1.0 + refreshDutyCycle);
	




            		 	    // colDelay = param->chargeDelay;
            		 	    // readLatency += colDelay;

            		 	    // // // --- DRAM Destructive Read Penalty (DRO) ---
            		 	    // // // Data is destroyed upon reading, requiring an immediate write-back pulse
            		 	    // // currentSenseAmp.CalculateLatency(columnResistance, 1, 1);
            		 	    // // double senseLatencyTotal = currentSenseAmp.readLatency * 1; 
            		 	    // // 
            		 	    // // // Find the slowest driver for the write-back phase
            		 	    // // double maxWriteDriverLatency = MAX(blSwitchMatrix.readLatency + blDecoder.readLatency, 
            		 	    // //                                    wlSwitchMatrix.readLatency + wlDecoder.readLatency);
            		 	    // // double writePhaseLatency = maxWriteDriverLatency + 10e-9; // 10ns write pulse width
            		 	    // // 
            		 	    // // // Add the Write-Back penalty to the total Read cycle
            		 	    // // readLatency += (senseLatencyTotal + writePhaseLatency);
            		 	    // // // -------------------------------------------

            		 	    // // double delay_BL = 0.69 * resCol * capBL;
            		 	    // // double delay_WL = 0.69 * resRow * capWL;
            		 	    // // readLatency += MAX(delay_BL, delay_WL);

            		 	    // // double blPathDelay = blDecoder.readLatency + blSwitchMatrix.readLatency;
            		 	    // // double wlPathDelay = wlDecoder.readLatency + wlSwitchMatrix.readLatency;
            		 	    // // readLatency += MAX(blPathDelay, wlPathDelay);

            		 	    // // readLatency += ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);
            		 	    // // readLatency += (param->SARADC ? sarADC.readLatency : (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency));
            		 	    // // readLatency *= (validated==true? param->beta : 1);

				    // double writeBackLatency = Trc * (numRow * activityRowRead);

        			    // // 1T1C DRAM Refresh Latency Penalty (Statistical Availability Model)
        			    // 
        			    // // 1. Retention Parameters
        			    // double tRefresh = 32e-3; // 32 ms retention time
        			    // double tRefreshPerRow = tRefresh / numRow; // 62.5 us for 512 rows

        			    // // 2. Latency of a Single Refresh Event
        			    // // A refresh only exercises the Array and Routing Drivers, NOT the ADCs or Accumulators.
        			    // // It consists of the peripheral delays + wire RC delays + the physical write-back pulse.
        			    // double refreshWriteBackPulse = 10e-9; // ~10ns physical charge restoration time
        			    // 
        			    // double refreshEventLatency = wlDecoder.readLatency + wlSwitchMatrix.readLatency 
        			    //                            + blDecoder.readLatency + blSwitchMatrix.readLatency 
        			    //                            + delay_setup // Wire RC charge-sharing delay
        			    //                            + refreshWriteBackPulse; 

        			    // // 3. Calculate Refresh Duty Cycle (Hardware Lockout Percentage)
        			    // // Example: If a refresh takes 20ns, and happens every 62,500ns, the array is locked 0.032% of the time.
        			    // double refreshDutyCycle = refreshEventLatency / tRefreshPerRow;

        			    // // 4. Apply Statistical Penalty to Total Read Latency
        			    // // By multiplying by (1 + DutyCycle), we statistically distribute the refresh collision stall-time 
        			    // // across all inference operations.
        			    // readLatency = readLatency * (1.0 + refreshDutyCycle);        
				    // 
				    // readLatency += writeBackLatency;
            		 	}
        		 }

        		 if (!CalculateclkFreq) {
        		     double maxDriverLatency = MAX(blSwitchMatrix.writeLatency + blDecoder.writeLatency, wlSwitchMatrix.writeLatency + wlDecoder.writeLatency);
        		     // writeLatency = (maxDriverLatency + 10e-9) * numWriteOperationPerRow * numRow;
			     
			     double writePulseTime = 10e-9;  
			    	/* NLS kinetics: tau(E) = tau_inf * exp[(Ea/E)^alpha].
				 * tau_inf = 236 ns is a HARD floor -- no field is faster.
				 * The old 10 ns was 24x below it.                       [4] */
			//	double eWrField = (cell.writeVoltage / param->tFE) / 1e8;
			//	double writePulseTime = param->nlsTauInf
			//	        * exp(pow(param->nlsEa / eWrField, param->nlsAlpha));
			     double numRowWritten = numRow * activityRowWrite;
                                writeLatency     = (maxDriverLatency + writePulseTime) * numWriteOperationPerRow * numRowWritten;
                                writeLatencyCore =  writePulseTime                     * numWriteOperationPerRow * numRowWritten;


        		     if (param->synchronous) {
        		         readLatencyADC = numColMuxed;
        		         readLatencyAccum = 0;
        		         if (numAdd > 1) readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);
        		     } else {
        		         double sensingDelay = (param->SARADC ? sarADC.readLatency : (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency));
        		         double driverDelay = MAX(blDecoder.readLatency + blSwitchMatrix.readLatency, wlDecoder.readLatency + wlSwitchMatrix.readLatency);
        		         driverDelay += ((numColMuxed > 1) ? (mux.readLatency + muxDecoder.readLatency) : 0);

        		         readLatencyADC = (driverDelay + sensingDelay) * numColMuxed * (validated==true? param->beta : 1);
        		         readLatencyAccum = MAX(((numAdd>1) ? adder.readLatency * numColMuxed*(numAdd-1) : 0) + ((numCellPerSynapse > 1) ? shiftAddWeight.readLatency : 0) + ((numReadPulse > 1) ? shiftAddInput.readLatency : 0) - readLatencyADC, 0.0);
        		     }

        		     if (numAdd > 1) {
        		         adder.CalculateLatency(1e20, dff.capTgDrain, 1);
        		         dff.CalculateLatency(1e20, 1);
        		     }
        		     if (numCellPerSynapse > 1) shiftAddWeight.CalculateLatency(numColMuxed);
        		     if (numReadPulse > 1) shiftAddInput.CalculateLatency(ceil(numColMuxed/numCellPerSynapse));

        		     readLatency = readLatencyADC + readLatencyAccum;
        		}

			 readLatency = 0;
			 readLatency = 60e-9;
		} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (conventionalSequential) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				double colRamp = 0;

				// 1.4 update: needs check
				double tau = (capCol)*(cell.resMemCellAvg);
				colDelay = horowitz(tau, 0, 1e20, &colRamp);	// Just to generate colRamp
				colDelay = tau * 0.2;  // assume the 15~20% voltage drop is enough for sensing
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {		
					
					// 1.4 update
					wlDecoder.CalculateLatency(1e20, capRow2, NULL, resRow, numCol, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					if (cell.accessType == CMOS_access) {
						wlNewDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow2, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);	
					} else {
						wlDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow1, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);										
					}					
					if (numColMuxed > 1) {
						mux.CalculateLatency(1e20, 0, 1);
						// 1.4 update
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						if (avgWeightBit > 1) {
							multilevelSAEncoder.CalculateLatency(1e20, 1);
						}
					}					
					if (CalculateclkFreq) {
						readLatency += MAX(wlDecoder.readLatency + wlNewDecoderDriver.readLatency + wlDecoderDriver.readLatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);		// latency factor of sensing cycle, beta = 1.4 by default					
					}
				}
				if (!CalculateclkFreq) {
					// Anni update
					adder.CalculateLatency(1e20, dff.capTgDrain, 1);	// numRead = numColMuxed*(numRow*activityRowRead-1)
					dff.CalculateLatency(1e20, 1);
					if (numCellPerSynapse > 1) {
						shiftAddWeight.CalculateLatency(1);	// numRead = (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse)
					}
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(1);	// numRead = ceil(numColMuxed/numCellPerSynapse)
					}
					if (param->synchronous) {
						readLatencyADC = numRow*activityRowRead*numColMuxed;
						// Anni update: hide readLatencyAccum by pipeline
						// adder is pipelined with ADC
						readLatencyAccum += numColMuxed*(numRow*activityRowRead-1) * (ceil(adder.readLatency*clkFreq)-1);
						if (numCellPerSynapse > 1) {	// shiftAddWeight is pipelined with adder+ADC of one column
							readLatencyAccum += (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddWeight.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/numColMuxed, 0);	
						} 
						if (numReadPulse > 1) {	 // shiftAddInput is pipelined with adder+ADC+shiftaddweight of numCellPerSynapse columns
							readLatencyAccum += ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddInput.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/ceil(numColMuxed/numCellPerSynapse), 0);	
						} 
					} else {
						readLatencyADC = (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency + colDelay) * (numRow*activityRowRead*numColMuxed) * (validated==true? param->beta : 1);
						readLatencyOther = MAX((wlDecoder.readLatency + wlNewDecoderDriver.readLatency + wlDecoderDriver.readLatency)*numRow*activityRowRead, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0)) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed*(numRow*activityRowRead-1) + shiftAddWeight.adder.readLatency * (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) + \
											shiftAddInput.adder.readLatency * ceil(numColMuxed/numCellPerSynapse) - readLatencyADC - readLatencyOther, 0);

					}
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
					// // Write
					// wllevelshifter.CalculateLatency(1e20, 2*wlNewDecoderDriver.capTgDrain, wlNewDecoderDriver.resTg, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					// bllevelshifter.CalculateLatency(1e20, 2*wlNewDecoderDriver.capTgDrain, wlNewDecoderDriver.resTg, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					// sllevelshifter.CalculateLatency(1e20, 2*slSwitchMatrix.capTgDrain, slSwitchMatrix.resTg, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					// slSwitchMatrix.CalculateLatency(1e20, capCol, resCol, 0, 2*numWriteOperationPerRow*numRow*activityRowWrite); 
					// writeLatencyArray = numWritePulse * param->writePulseWidthLTP + (-numErasePulse) * param->writePulseWidthLTD;
					// writeLatency += MAX(wlDecoder.writeLatency + wlNewDecoderDriver.writeLatency + wlDecoderDriver.writeLatency, sllevelshifter.writeLatency + slSwitchMatrix.writeLatency + bllevelshifter.writeLatency);
					// writeLatency += writeLatencyArray;
					
			} else if (conventionalParallel) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double colRamp = 0;

				// 1.4 update: needs check - (capCol)*(cell.resMemCellAvg/(numRow));? 	Anni update: numRow->numRowParallel
				double tau = (capCol)*(cell.resMemCellAvg/(numRowParallel/2));

				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2;  // assume the 15~20% voltage drop is enough for sensing
				if (CalculateclkFreq || !param->synchronous) {				
					double bufferlatency=0;			 
					if (cell.accessType == CMOS_access) {
						// 1.4 update : buffer latency
						int iterbuffnum = param->buffernumber-1;
						if (param->buffernumber>0) {
							double buffnum = param->buffernumber; // # of buffers
							double sectionnum = param->numColSubArray/(param->buffernumber+1); // # of cells in each section
							double unitcap= capRow2/param->numColSubArray;
							double unitres= resRow/param->numColSubArray;

							param->unitcap = unitcap;
							param->unitres = unitres;
							param->drivecapin = drivecapin;

							double capload_repeater1 = capRow2/(buffnum+1)+ drivecapin;
							double capload_repeater2 = capRow2/(buffnum+1);
							
							wlNewSwitchMatrix.CalculateLatency(1e20, capload_repeater1, unitres * sectionnum, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
							
							double t1 =  targetdriveres * (drivecapout) * 0.69 
										+ unitcap * sectionnum * (0.69*targetdriveres + 0.38* unitres * sectionnum)
										+ (unitcap * sectionnum * 0.69 + 0.69 * targetdriveres)* drivecapin;											
							double t2 =  targetdriveres * (drivecapout) * 0.69 
										+ unitcap * sectionnum * (0.69*targetdriveres + 0.38* unitres * sectionnum);									
							double t3 = (drivecapout + drivecapin) * targetdriveres * 0.69;				
							while (iterbuffnum >= 0){
								if (iterbuffnum  == 0) {
									bufferlatency += t2 + t3;
								} else {
									bufferlatency += t1 + t3;
								}
								iterbuffnum  = iterbuffnum-1;
							}
						} else {
							bufferlatency=0;
							wlNewSwitchMatrix.CalculateLatency(1e20, capRow2, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
						}
					} else {
						wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					}
					if (numColMuxed>1) {
						mux.CalculateLatency(1e20, 0, 1);
						// 1.4 update
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}				
					if (CalculateclkFreq) {
						// 1.4 update - updated
						readLatency += MAX(wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency + bufferlatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						// readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default									
						param->rowdelay = wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency + bufferlatency;
						param->muxdelay = mux.readLatency+muxDecoder.readLatency;
						param->ADClatency = multilevelSenseAmp.readLatency;					
					
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: add partial sums
					if (numAdd > 1) {	
						adder.CalculateLatency(1e20, dff.capTgDrain, 1); // numRead = numColMuxed*(numAdd-1)
						dff.CalculateLatency(1e20, 1);	// numRead = numColMuxed*(numAdd-1)
					}
					if (numCellPerSynapse > 1) {
						shiftAddWeight.CalculateLatency(1);	// numRead = (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse)
					}
					if (numReadPulse > 1) {
						shiftAddInput.CalculateLatency(1);	// numRead = ceil(numColMuxed/numCellPerSynapse)
					}
					if (param->synchronous) {
						readLatencyADC = numColMuxed * numAdd;	// Anni update	
						// Anni update: hide readLatencyAccum by pipeline
						if (numAdd > 1) {	// adder is pipelined with ADC
							readLatencyAccum += numColMuxed*(numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);	
						}	
						if (numCellPerSynapse > 1) {	// shiftAddWeight is pipelined with adder+ADC of one column
							readLatencyAccum += (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddWeight.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/numColMuxed, 0);	
						} 
						if (numReadPulse > 1) {	 // shiftAddInput is pipelined with adder+ADC+shiftaddweight of numCellPerSynapse columns
							readLatencyAccum += ceil(numColMuxed/numCellPerSynapse) * MAX(ceil(shiftAddInput.adder.readLatency*clkFreq) - (readLatencyADC+readLatencyAccum)/ceil(numColMuxed/numCellPerSynapse), 0);	
						} 
					} else {
						readLatencyADC = (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency + colDelay) * numColMuxed * (validated==true? param->beta : 1) * numAdd;
						readLatencyOther = MAX((wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency) * numAdd, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0)) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed*(numAdd-1) + shiftAddWeight.adder.readLatency * (numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse) + \
											shiftAddInput.adder.readLatency * ceil(numColMuxed/numCellPerSynapse) - readLatencyADC - readLatencyOther, 0);
					}
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}

				 // cout << "\n================ RRAM/FeFET Block Latency ================" << endl;
        			 // cout << "--- Read Path ---" << endl;
        			 // if (wlDecoder.initialized) cout << "WL Decoder Read Latency:       " << wlDecoder.readLatency << " ns" << endl;
        			 // if (wlNewDecoderDriver.initialized) cout << "WL New Decoder Driver Read Latency: " << wlNewDecoderDriver.readLatency << " ns" << endl;
        			 // if (wlDecoderDriver.initialized) cout << "WL Decoder Driver Read Latency:" << wlDecoderDriver.readLatency << " ns" << endl;
        			 // if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Read Latency: " << wlSwitchMatrix.readLatency << " ns" << endl;
        			 // if (wlNewSwitchMatrix.initialized) cout << "WL New Switch Matrix Read Latency: " << wlNewSwitchMatrix.readLatency << " ns" << endl;
        			 // if (slSwitchMatrix.initialized) cout << "SL Switch Matrix Read Latency: " << slSwitchMatrix.readLatency << " ns" << endl;

        			 // if (numColMuxed > 1) {
        			 //     if (muxDecoder.initialized) cout << "MUX Decoder Read Latency:      " << muxDecoder.readLatency << " ns" << endl;
        			 //     if (mux.initialized) cout << "MUX Read Latency:              " << mux.readLatency << " ns" << endl;
        			 // }
        			 // if (param->SARADC) {
        			 //     if (sarADC.initialized) cout << "SAR ADC Read Latency:          " << sarADC.readLatency << " ns" << endl;
        			 // } else {
        			 //     if (multilevelSenseAmp.initialized) cout << "Multilevel SA Read Latency:    " << multilevelSenseAmp.readLatency << " ns" << endl;
        			 //     if (multilevelSAEncoder.initialized) cout << "Multilevel SA Encoder Latency: " << multilevelSAEncoder.readLatency << " ns" << endl;
        			 // }
        			 // if (numAdd > 1) {
        			 //     if (adder.initialized) cout << "Adder Read Latency:            " << adder.readLatency << " ns" << endl;
        			 //     if (dff.initialized) cout << "DFF Read Latency:              " << dff.readLatency << " ns" << endl;
        			 // }
        			 // if (numCellPerSynapse > 1) {
        			 //     if (shiftAddWeight.initialized) cout << "ShiftAdd Weight Read Latency:  " << shiftAddWeight.readLatency << " ns" << endl;
        			 // }
        			 // if (numReadPulse > 1) {
        			 //     if (shiftAddInput.initialized) cout << "ShiftAdd Input Read Latency:   " << shiftAddInput.readLatency << " ns" << endl;
        			 // }

        			 // cout << "--- Write Path ---" << endl;
        			 // if (wlDecoder.initialized) cout << "WL Decoder Write Latency:      " << wlDecoder.writeLatency << " ns" << endl;
        			 // if (wlNewDecoderDriver.initialized) cout << "WL New Decoder Driver Write Latency: " << wlNewDecoderDriver.writeLatency << " ns" << endl;
        			 // if (wlDecoderDriver.initialized) cout << "WL Decoder Driver Write Latency:" << wlDecoderDriver.writeLatency << " ns" << endl;
        			 // if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Write Latency:" << wlSwitchMatrix.writeLatency << " ns" << endl;
        			 // if (wlNewSwitchMatrix.initialized) cout << "WL New Switch Matrix Write Latency:" << wlNewSwitchMatrix.writeLatency << " ns" << endl;
        			 // if (slSwitchMatrix.initialized) cout << "SL Switch Matrix Write Latency:" << slSwitchMatrix.writeLatency << " ns" << endl;

        			 // if (cell.writeVoltage > 1.5) {
        			 //     if (wllevelshifter.initialized) cout << "WL Level Shifter Write Latency:" << wllevelshifter.writeLatency << " ns" << endl;
        			 //     if (bllevelshifter.initialized) cout << "BL Level Shifter Write Latency:" << bllevelshifter.writeLatency << " ns" << endl;
        			 //     if (sllevelshifter.initialized) cout << "SL Level Shifter Write Latency:" << sllevelshifter.writeLatency << " ns" << endl;
        			 // }
        			 // cout << "==========================================================" << endl;
			} else if (BNNsequentialMode || XNORsequentialMode) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				double colRamp = 0;

				// 1.4 update: needs check
				double tau = (capCol)*(cell.resMemCellAvg);
				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2 ;  // assume the 15~20% voltage drop is enough for sensing
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				if (CalculateclkFreq || !param->synchronous) {

					// 1.4 update
					wlDecoder.CalculateLatency(1e20, capRow2, NULL, resRow, numCol, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					if (cell.accessType == CMOS_access) {
						wlNewDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow2, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);	
					} else {
						wlDecoderDriver.CalculateLatency(wlDecoder.rampOutput, capRow1, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					}
					if (numColMuxed > 1) {
						mux.CalculateLatency(1e20, 0, 1);
						// 1.4 update
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}
						// 1.4 update 230615
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);

					if (CalculateclkFreq) {
						readLatency += MAX(wlDecoder.readLatency + wlNewDecoderDriver.readLatency + wlDecoderDriver.readLatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						// readLatency += colDelay;
						// 1.4 update 230615
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				}
				if (!CalculateclkFreq) {
					// Anni update: hide readLatencyAccum by pipeline
					adder.CalculateLatency(1e20, dff.capTgDrain, 1);
					dff.CalculateLatency(1e20, 1);
					if (param->synchronous) {
						readLatencyADC = numRow*activityRowRead*numColMuxed;
						// adder is pipelined with ADC
						readLatencyAccum = numColMuxed*(numRow*activityRowRead-1) * (ceil(adder.readLatency*clkFreq)-1);		
					} else { 
						readLatencyADC = (rowCurrentSenseAmp.readLatency + colDelay) * numRow*activityRowRead*numColMuxed * (validated==true? param->beta : 1);
						readLatencyOther = MAX((wlDecoder.readLatency + wlNewDecoderDriver.readLatency + wlDecoderDriver.readLatency)*numRow*activityRowRead, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0)) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed*(numRow*activityRowRead-1) - readLatencyADC - readLatencyOther, 0);
					}					
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
			} else if (BNNparallelMode || XNORparallelMode) {
				double capBL = lengthCol * 0.2e-15/1e-6;
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double colRamp = 0;

				// 1.4 update: needs check; Anni update: numRow->numRowParallel
				double tau = (capCol)*(cell.resMemCellAvg/(numRowParallel/2));
				colDelay = horowitz(tau, 0, 1e20, &colRamp);
				colDelay = tau * 0.2;  // assume the 15~20% voltage drop is enough for sensing
				if (CalculateclkFreq || !param->synchronous) {
					if (cell.accessType == CMOS_access) {
						wlNewSwitchMatrix.CalculateLatency(1e20, capRow2, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					} else {
						wlSwitchMatrix.CalculateLatency(1e20, capRow1, resRow, 1, 2*numWriteOperationPerRow*numRow*activityRowWrite);
					}
					if (numColMuxed > 1) {
						mux.CalculateLatency(1e20, 0, 1);
						// 1.4 update 
						muxDecoder.CalculateLatency(1e20, mux.capTgGateN*ceil(numCol/numColMuxed), mux.capTgGateP*ceil(numCol/numColMuxed), 0, 0, 1, 0);
					}
					if (param->SARADC) {
						sarADC.CalculateLatency(1);
					} else {
						multilevelSenseAmp.CalculateLatency(columnResistance, 1, 1);
						multilevelSAEncoder.CalculateLatency(1e20, 1);
					}
					if (CalculateclkFreq) {
						readLatency += MAX(wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0));
						// readLatency += colDelay;
						readLatency += multilevelSenseAmp.readLatency;
						readLatency += multilevelSAEncoder.readLatency;
						readLatency += sarADC.readLatency;
						readLatency *= (validated==true? param->beta : 1);	// latency factor of sensing cycle, beta = 1.4 by default
					}
				}
				if (!CalculateclkFreq) {					
					// Anni update: add partial sums; hide readLatencyAccum by pipeline
					if (numAdd > 1) {	
						adder.CalculateLatency(1e20, dff.capTgDrain, 1);
						dff.CalculateLatency(1e20, 1);
					}
					if (param->synchronous) {
						readLatencyADC = numColMuxed * numAdd;
						// adder is pipelined with ADC
						readLatencyAccum = numColMuxed * (numAdd-1) * (ceil(adder.readLatency*clkFreq)-1);
					} else { 
						readLatencyADC = (multilevelSenseAmp.readLatency + multilevelSAEncoder.readLatency + sarADC.readLatency + colDelay) * numColMuxed * (validated==true? param->beta : 1) * numAdd;
						readLatencyOther = MAX((wlNewSwitchMatrix.readLatency + wlSwitchMatrix.readLatency) * numAdd, ((numColMuxed > 1)==true? (mux.readLatency+muxDecoder.readLatency):0)) * numColMuxed * (validated==true? param->beta : 1);
						// Anni update: hide readLatencyAccum by pipeline
						readLatencyAccum = MAX(adder.readLatency * numColMuxed * (numAdd-1) - readLatencyADC - readLatencyOther, 0);
					}					
					readLatency = readLatencyADC + readLatencyAccum + readLatencyOther;
				}
			}
		}
	}
}

void SubArray::CalculatePower(const vector<double> &columnResistance, bool writeBack) {
	readEnergyCore = readSenseEnergyCore = readRestoreEnergyCore = 0;
	writeEnergyCore = writeCellEnergyCore = inhibitionEnergyCore = 0;

	if (!initialized) {
		cout << "[Subarray] Error: Require initialization first!" << endl;
	} else {
		readDynamicEnergy = 0;
		writeDynamicEnergy = 0;
		readDynamicEnergyArray = 0;
		
		double numReadOperationPerRow;   // average value (can be non-integer for energy calculation)
		if (numCol > numReadCellPerOperationNeuro)
			numReadOperationPerRow = numCol / numReadCellPerOperationNeuro;
		else
			numReadOperationPerRow = 1;

		double numWriteOperationPerRow;   // average value (can be non-integer for energy calculation)
		if (numCol * activityColWrite > numWriteCellPerOperationNeuro)
			numWriteOperationPerRow = numCol * activityColWrite / numWriteCellPerOperationNeuro;
		else
			numWriteOperationPerRow = 1;

		if (cell.memCellType == Type::SRAM) {
			
			// Array leakage (assume 2 INV)
			leakage = 0;
			leakage += CalculateGateLeakage(INV, 1, cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize,
					cell.widthSRAMCellPMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, inputParameter.temperature, tech) * tech.vdd * 2;
			// Anni update
			leakageSRAMInUse = leakage;
			leakage *= numRow * numCol;

			if (conventionalSequential) {
				wlDecoder.CalculatePower(numRow*activityRowRead, numRow*activityRowWrite);
				precharger.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numWriteOperationPerRow*numRow*activityRowWrite);
				sramWriteDriver.CalculatePower(numWriteOperationPerRow*numRow*activityRowWrite);
				adder.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numReadCellPerOperationNeuro/numCellPerSynapse);				
				// Anni update
				dff.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numReadCellPerOperationNeuro/numCellPerSynapse*(ceil(log2(numRow))/2+1), param->validated);
				senseAmp.CalculatePower(numReadOperationPerRow*numRow*activityRowRead);
				// Anni update:
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculatePower(numCellPerSynapse-1);	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(1);					
				}

				// Array
				// 1.4 update: read energy update
				readDynamicEnergyArray = capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowRead;  // Just BL discharging // -added, wordline charging
				
				// // 1.4 update: WL energy for write + modification (assuming toggling of SRAM bit at every write, each Q/Qbar consumes half CVdd^2)
				// writeDynamicEnergyArray = cell.capSRAMCell * tech.vdd * tech.vdd * numCol * activityColWrite * numRow * activityRowWrite;    // flip Q and Q_bar
				// writeDynamicEnergyArray += capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowWrite;
				
				// Read
				readDynamicEnergy += wlDecoder.readDynamicEnergy;
				readDynamicEnergy += precharger.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += senseAmp.readDynamicEnergy;
				readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				
				readDynamicEnergyADC = precharger.readDynamicEnergy + readDynamicEnergyArray + senseAmp.readDynamicEnergy;
				readDynamicEnergyAccum = adder.readDynamicEnergy + dff.readDynamicEnergy + shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergyOther = wlDecoder.readDynamicEnergy;

				// Write
				// writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
				// writeDynamicEnergy += precharger.writeDynamicEnergy;
				// writeDynamicEnergy += sramWriteDriver.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;
				
				// Leakage
				leakage += wlDecoder.leakage;
				leakage += precharger.leakage;
				leakage += sramWriteDriver.leakage;
				leakage += senseAmp.leakage;
				leakage += dff.leakage;
				leakage += adder.leakage;
				leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
				// Anni update
				leakageSRAMInUse *= (numRow-1) * numCol;

			} else if (conventionalParallel) {
				wlSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				// Anni update
				precharger.CalculatePower(numColMuxed*numAdd, numWriteOperationPerRow*numRow*activityRowWrite);
				sramWriteDriver.CalculatePower(numWriteOperationPerRow*numRow*activityRowWrite);
				
				// 1.4 update: ADC update
				param->reference_energy_peri = capRow1/param->numColSubArray * tech.vdd * tech.vdd * (numRow);				
				
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				// Anni update: numAdd
				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numAdd);
				} else {
					multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
					multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
				}
				if (numAdd > 1) {	
					adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)); 
					dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculatePower((numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse));	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
				}
				// Array

				// 1.4 update: read energy update; Anni update: * numColMuxed
				readDynamicEnergyArray = capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowRead; // added for WL/BL discharging
				// 1.4 update: buffer energy
				readDynamicEnergyArray += (drivecapin + drivecapout) * tech.vdd * tech.vdd * param->buffernumber * 2 * numRow * activityRowRead;
				// 1.4 update: iterate for numColMuxed - needs check if this is necessary
				readDynamicEnergyArray *= numColMuxed;
				
				// // 1.4 update: WL energy for write + modification (assuming toggling of SRAM bit at every write, each Q/Qbar consumes half CVdd^2)
				// writeDynamicEnergyArray = cell.capSRAMCell * tech.vdd * tech.vdd * numCol * activityColWrite * numRow * activityRowWrite;    // flip Q and Q_bar
				// writeDynamicEnergyArray += capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowWrite;

				
				// Read
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				// readDynamicEnergy += precharger.readDynamicEnergy; -> precharger is not needed for SRAM parallel mode
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				// Anni update
				readDynamicEnergy += ((numColMuxed > 1)==true? mux.readDynamicEnergy:0);
				readDynamicEnergy += ((numColMuxed > 1)==true? muxDecoder.readDynamicEnergy:0);
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergy += sarADC.readDynamicEnergy;

				// 1.4 update : precharger not needed 230619
				readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;				
				// Anni update
				readDynamicEnergyAccum = shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy + adder.readDynamicEnergy + dff.readDynamicEnergy;
				readDynamicEnergyOther = wlSwitchMatrix.readDynamicEnergy + ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				
				// Write
				// writeDynamicEnergy += wlSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += precharger.writeDynamicEnergy;
				// writeDynamicEnergy += sramWriteDriver.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage += wlSwitchMatrix.leakage;
				leakage += precharger.leakage;
				leakage += sramWriteDriver.leakage;
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
				// Anni update
				leakage += dff.leakage;
				leakage += adder.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				// Anni update
				leakageSRAMInUse *= (numRow * numCol - numRowParallel * ceil(numCol/numColMuxed));
			
			} else if (BNNsequentialMode || XNORsequentialMode) {
				wlDecoder.CalculatePower(numRow*activityRowRead, numRow*activityRowWrite);
				precharger.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numWriteOperationPerRow*numRow*activityRowWrite);
				sramWriteDriver.CalculatePower(numWriteOperationPerRow*numRow*activityRowWrite);
				// Anni update: 
				adder.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numReadCellPerOperationNeuro);				
				dff.CalculatePower(numReadOperationPerRow*numRow*activityRowRead, numReadCellPerOperationNeuro*(ceil(log2(numRow))/2+1), param->validated);
				senseAmp.CalculatePower(numReadOperationPerRow*numRow*activityRowRead);
				
				// 1.4 update: read energy update : needs check for BNN/XNOR mode
				readDynamicEnergyArray = capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowRead; // added for WL/BL discharging
				
				// // 1.4 update: WL energy for write + modification (assuming toggling of SRAM bit at every write, each Q/Qbar consumes half CVdd^2)
				// // needs check for BNN/XNOR mode
				// writeDynamicEnergyArray = cell.capSRAMCell * tech.vdd * tech.vdd * numCol * activityColWrite * numRow * activityRowWrite;    // flip Q and Q_bar
				// writeDynamicEnergyArray += capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowWrite;

				// Read
				readDynamicEnergy += wlDecoder.readDynamicEnergy;
				readDynamicEnergy += precharger.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += senseAmp.readDynamicEnergy;
				
				// Write				
				// writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
				// writeDynamicEnergy += precharger.writeDynamicEnergy;
				// writeDynamicEnergy += sramWriteDriver.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage += wlDecoder.leakage;
				leakage += precharger.leakage;
				leakage += sramWriteDriver.leakage;
				leakage += senseAmp.leakage;
				leakage += dff.leakage;
				leakage += adder.leakage;
				// Anni update
				leakageSRAMInUse *= (numRow-1) * numCol;
				
			} else if (BNNparallelMode || XNORparallelMode) {
				wlSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				// Anni update
				precharger.CalculatePower(numColMuxed * numAdd, numWriteOperationPerRow*numRow*activityRowWrite);
				sramWriteDriver.CalculatePower(numWriteOperationPerRow*numRow*activityRowWrite);
				// Anni update: add mux; numAdd
				if (numColMuxed>1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numAdd);
				} else {
					multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
					multilevelSAEncoder.CalculatePower(numColMuxed*numAdd);
				}
				if (numAdd > 1) {	
					adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)); 
					dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
				}
				// Array
				// 1.4 update: read energy update : needs check for BNN/XNOR mode; Anni update
				readDynamicEnergyArray = capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowRead * numColMuxed; // added for WL/BL discharging
				// 1.4 update: ADC update
				param->reference_energy_peri = capRow1/param->numColSubArray * tech.vdd * tech.vdd * numRow;
				
				// // 1.4 update: WL energy for write + modification (assuming toggling of SRAM bit at every write, each Q/Qbar consumes half CVdd^2)
				// // needs check for BNN/XNOR mode
				// writeDynamicEnergyArray = cell.capSRAMCell * tech.vdd * tech.vdd * numCol * activityColWrite * numRow * activityRowWrite;    // flip Q and Q_bar
				// writeDynamicEnergyArray += capRow1 * tech.vdd * tech.vdd * (numRow) * activityRowWrite;

				// Read
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += precharger.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += sarADC.readDynamicEnergy;
				// Anni update:
				readDynamicEnergy += ((numColMuxed > 1)==true? mux.readDynamicEnergy:0);
				readDynamicEnergy += ((numColMuxed > 1)==true? muxDecoder.readDynamicEnergy:0);
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				
				// Write				
				// writeDynamicEnergy += wlSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += precharger.writeDynamicEnergy;
				// writeDynamicEnergy += sramWriteDriver.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage += wlSwitchMatrix.leakage;
				leakage += precharger.leakage;
				leakage += sramWriteDriver.leakage;
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				// Anni update
				leakage += dff.leakage;
				leakage += adder.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				// Anni update
				leakageSRAMInUse *= (numRow * numCol - numRowParallel * ceil(numCol/numColMuxed));
				
			}		
	    } 
		else if (cell.memCellType == Type::Cap) { //nvCap added
			leakageSRAMInUse = 0;
			if (conventionalParallel) {
				double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double capBL = lengthCol * 0.2e-15/1e-6;
			
				wlSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				// Anni update: numAdd
				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numAdd);
				} else {
					// multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
					multilevelSenseAmp.readDynamicEnergy = 0;
					double P_Col = 1.7*(levelOutput-1)*1e-6/2;
					multilevelSenseAmp.readDynamicEnergy += P_Col*1e-9 * columnResistance.size();
					multilevelSenseAmp.readDynamicEnergy *= numAdd;
					multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
				}
				if (numAdd > 1) {	
					adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)); 
					dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculatePower((numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse));	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
				}
				// Read
				readDynamicEnergyArray = 0;
				// Anni update:  * numAdd
				// readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells * numAdd; // Selected BLs activityColWrite
				readDynamicEnergyArray += capRow2 * cell.readVoltage * cell.readVoltage * numRow * activityRowRead; // Selected WL
				readDynamicEnergyArray *= numColMuxed;
				readDynamicEnergyArray += 0.29e-3 * (colDelay + 1e-9) * numCol / 128 * numAdd;

				readDynamicEnergy = 0;
				readDynamicEnergy += wlNewSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				// Anni update: adder, dff, mux
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += sarADC.readDynamicEnergy;	
				
				readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;
				// Anni update: accum, other
				readDynamicEnergyAccum = shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy + adder.readDynamicEnergy + dff.readDynamicEnergy;				
				readDynamicEnergyOther = wlNewSwitchMatrix.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy + ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);								
				
				// Leakage
				leakage = 0;
				leakage += wlSwitchMatrix.leakage;
				leakage += wlNewSwitchMatrix.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
				// Anni update
				leakage += dff.leakage;
				leakage += adder.leakage;

				// cout<<"wlSwitchMatrix.leakage: "<<wlSwitchMatrix.leakage<<endl;
				// cout<<"wlNewSwitchMatrix.leakage: "<<wlNewSwitchMatrix.leakage<<endl;
				// cout<<"slSwitchMatrix.leakage: "<<slSwitchMatrix.leakage<<endl;
				// cout<<"mux.leakage: "<<mux.leakage<<endl;
				// cout<<"muxDecoder.leakage: "<<muxDecoder.leakage<<endl;
				// cout<<"multilevelSenseAmp.leakage: "<<multilevelSenseAmp.leakage<<endl;
				// cout<<"multilevelSAEncoder.leakage: "<<multilevelSAEncoder.leakage<<endl;
				// cout<<"shiftAddWeight.leakage: "<<shiftAddWeight.leakage<<endl;
				// cout<<"shiftAddInput.leakage: "<<shiftAddInput.leakage<<endl;
				// cout<<"dff.leakage: "<<dff.leakage<<endl;
				// cout<<"adder.leakage: "<<adder.leakage<<endl;
				// cout<<"leakage: "<<leakage<<endl;
			}
		

	
		 } else if (cell.memCellType == Type::_2TnC) {
                	double bitsPerCell = param->bitsPerCell;
			double numReadOperationPerRow = (numCol > numReadCellPerOperationNeuro) ? numCol / numReadCellPerOperationNeuro : 1;
			double totalCellFlippingEnergy = 0;

                	//Calculate Operation Counts
                	//How many times it loops to process a single logical row
                	double numReadCells = (int)ceil((double)numCol/numColMuxed);
                	int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);

                	//Calculate Driver Power (Switching Energy of wires)

		
			// allel write requires the peripheral drivers to fire exactly twice per cell write
			//int writeCyclesPerCell = 2;
			double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite;


                	//RSL: Active during Read
                	rslDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, numRow*activityRowWrite);
                	rslSwitchMatrix.CalculatePower(numColMuxed, 0, activityRowRead, activityColWrite);

                	//WBL (Write Bitline - Row Plane): Active during Write
                	//Drives 8 wires, but activity is based on logical row access
			// wblDecoder.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite);
                	// wblSwitchMatrix.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);

			// wblDecoder.CalculatePower(0, totalWriteActivations);
                        // wblSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

			double wblReadActivations = numRow * activityRowRead * numColMuxed;

			wblDecoder.CalculatePower(wblReadActivations, totalWriteActivations);
			wblPlaneDecoder.CalculatePower(wblReadActivations, totalWriteActivations);
			// both decoders fire on every access, so energies ADD
			// (latencies are parallel and take MAX -- see the latency edit)
			wblDecoder.readDynamicEnergy  += wblPlaneDecoder.readDynamicEnergy;
			wblDecoder.writeDynamicEnergy += wblPlaneDecoder.writeDynamicEnergy;
			wblDecoder.leakage            += wblPlaneDecoder.leakage;
			wblSwitchMatrix.CalculatePower(numColMuxed, totalWriteActivations, activityRowRead, activityColWrite);


                	//SSL: Active during Write
                	// sslDecoder.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite);
                	// sslSwitchMatrix.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);

			sslDecoder.CalculatePower(0, totalWriteActivations);
                        sslSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);


                	//WWL: Active during Write
                	// wwlDecoder.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite);
                	// wwlSwitchMatrix.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);

			wwlDecoder.CalculatePower(0, totalWriteActivations);
                        wwlSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                	//WPL: Active during Write
                	// wplDecoder.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite);
                	// wplSwitchMatrix.CalculatePower(0, numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);

			wplDecoder.CalculatePower(0, totalWriteActivations);
                        wplSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                	//RBL: Via Mux
                	rblDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, 0);
                	rblSwitchMatrix.CalculatePower(numColMuxed, 0, activityRowRead, activityColWrite);


			// 	// 2-Cycle parallel write requires the peripheral drivers to fire exactly twice per cell write
		// 	// int writeCyclesPerCell = 2;
		// 	double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite;


                // 	//RSL (Col): Active during Read
		// 	rslDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, 0);
		// 	rslSwitchMatrix.CalculatePower(numRow*activityRowRead*numColMuxed, 0, activityRowRead, activityColWrite);

                // 	//WBL (Write Bitline - Row Plane): Active during Write
		// 	wblDecoder.CalculatePower(0, totalWriteActivations);
		// 	wblSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);
		// 	
		// 	//SSL (Col): Active during Write
		// 	sslDecoder.CalculatePower(0, totalWriteActivations);
		// 	sslSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                // 	//WWL (Row): Active during Write
		// 	wwlDecoder.CalculatePower(0, totalWriteActivations);
                //         wwlSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                // 	//WPL (Row): Active during Write
		// 	wplDecoder.CalculatePower(0, totalWriteActivations);
                //         wplSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                // 	// Calculate total number of rows being actively read
		// 	double totalReadActivations = numRow * activityRowRead;
		// 	
		// 	//RBL (Row): Via Mux
		// 	rblDecoder.CalculatePower(totalReadActivations * numColMuxed, 0);
		// 	rblSwitchMatrix.CalculatePower(totalReadActivations * numColMuxed, 0, activityRowRead, activityColWrite);
			
			//PERIPHERALS (Mux, ADC, etc)
			if (numColMuxed > 1) {
				mux.CalculatePower(numColMuxed);
				muxDecoder.CalculatePower(numColMuxed, 1);
			}
			
			if (param->SARADC) {
				sarADC.CalculatePower(columnResistance, 1);
			} else {
				multilevelSenseAmp.CalculatePower(columnResistance, 1);
				multilevelSAEncoder.CalculatePower(numColMuxed);
			}
			
			//Accumulation Logic
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculatePower(numColMuxed);
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));
			}
			if (numAdd > 1) {
				adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed));
				dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
			}
			

                	//Calculate Read Dynamic Energy (Array + Drivers)
			

			// Map 3D line parasitics to NeuroSim's calculated array parasitics
			double resCellWBL = resRow;
			
			double resCellRBL = resRow;
			
			double resCellRSL = resCol;
			
			double resReadTransistor = cell.resMemCellOn; // Assuming worst-case ON resistance

			double e_WBL_ON = 0;
                        double e_WBL_OFF = 0;
			double ReadEnergyArray = 0;
			double dielectric = 0;
			
			double e_RBL_ON = 0;
                        double e_RBL_OFF = 0;
                        double eReadOn = 0;
                        double eReadOff = 0;

			// double iFlip_ON_W = 20e-6;
                        // double tWrite = 10e-9;
                        // double CapFe = 0.5e-15;

                        // double eWritePerCell = 0;

						double iFlip_ON_W = 20e-6;
                        double tWrite = 10e-9;

                        double CapFe = EPS0_VACUUM * param->epsFE * param->cellAreaFE / param->tFE;

			double capGate = CapFe / param->capGateRatio;
			double capNode = param->bitsPerCell * CapFe + capGate;
			double QswRd   = param->twoPr * param->cellAreaFE;

			double vIntBias = param->vthReadTr + param->readOverdrive;
			double vUnsel   = (param->bitsPerCell > 1)
			                ? (vIntBias * capNode - CapFe * cell.readVoltage)
			                  / ((param->bitsPerCell - 1) * CapFe)
			                : 0.0;
			double dVsense  = QswRd / capNode;

			double iRead1 = MIN(param->ioffReadTr
			                    * pow(10.0, (vIntBias + 0.5*dVsense - param->vthReadTr) / param->ssReadTr),
			                    param->ionSatReadTr);
			double iRead0 = MIN(param->ioffReadTr
			                    * pow(10.0, (vIntBias - 0.5*dVsense - param->vthReadTr) / param->ssReadTr),
			                    param->ionSatReadTr);

			double p1Rd = param->dataOnesRead;

                        double eWritePerCell = 0;

			
			if (cell.mem_rdo == Type::ndro) {
				
				double alphaRd = 0.0;
				double V_RBL   = 0.5;
				double tREAD   = 10e-9;

				dielectric = CapFe * cell.readVoltage * cell.readVoltage;

				e_WBL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_WBL_OFF = 0.0;

				e_RBL_ON  = V_RBL * iRead1 * tREAD;
				e_RBL_OFF = V_RBL * iRead0 * tREAD;

				eReadOn  = e_WBL_ON  + e_RBL_ON  + dielectric;
				eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;
				

				// // 1. WBL Energy: Driven by the input drivers (DACs/PWM)
        			// // V_read is applied to the WBL of the target cells
        			// 
				// double iFlip_ON = 2e-6; // 2 uA instantaneous switching current
				// double iFlip_OFF = 1.5e-6; 
				// double tFlip = 0.7e-9; 
				// 
				// double V_RBL = 0.5; 
				// double iRBL_ON = 2e-6; 
				// double iRBL_OFF = 10e-9; 
				// double tREAD = 10e-9; 
				// 
				// dielectric = 0.5 * CapFe * cell.readVoltage * cell.readVoltage;
				// 
				// // 1. WBL:
				// e_WBL_ON = cell.readVoltage * iFlip_ON * tFlip;
        			// e_WBL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

        			// // 2. RBL Energy: Driven by a small bias voltage during read
				// e_RBL_ON = V_RBL * iRBL_ON * tREAD;
                        	// e_RBL_OFF = V_RBL * iRBL_OFF * tREAD;

				// eReadOn = e_WBL_ON + e_RBL_ON + dielectric;
				// eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

				// // ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow) + (0.5 * eReadOff * activityRowRead * numRow);

				// ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                   (0.5 * eReadOff * activityRowRead * numRow * numCol);
			
			} else if (cell.mem_rdo == Type::qndro) {

				double alphaRd = 0.20;
				double V_RBL   = 0.5;
				double tREAD   = 10e-9;

				dielectric = CapFe * cell.readVoltage * cell.readVoltage;

				e_WBL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_WBL_OFF = 0.0;

				e_RBL_ON  = V_RBL * iRead1 * tREAD;
				e_RBL_OFF = V_RBL * iRead0 * tREAD;

				eReadOn  = e_WBL_ON  + e_RBL_ON  + dielectric;
				eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;
				
				// // 1. WBL Energy: Driven by the input drivers (DACs/PWM)
                                // // V_read is applied to the WBL of the target cells

				// double iFlip_ON = 10e-6; // 20 uA instantaneous switching current
                                // double iFlip_OFF = 1.5e-6;
                                // double tFlip = 2e-9;

                                // double V_RBL = 0.5;
                                // double iRBL_ON = 2e-6; // NEED TO REDUCE
                                // double iRBL_OFF = 10e-9;
                                // double tREAD = 10e-9;

                                // dielectric = 0.5 * CapFe * cell.readVoltage * cell.readVoltage;
				// 
				// // 1. WBL Energy:
				// e_WBL_ON = cell.readVoltage * iFlip_ON * tFlip;
                                // e_WBL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

                                // // 2. RBL Energy: Driven by a small bias voltage during read
                                // e_RBL_ON = V_RBL * iRBL_ON * tREAD;
                                // e_RBL_OFF = V_RBL * iRBL_OFF * tREAD;

                                // eReadOn = e_WBL_ON + e_RBL_ON + dielectric;
                                // eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

                                // //double ReadEnergyArray = (0.5 * eReadOn * numReadCellPerOperationNeuro) +
                                // //                         (0.5 * eReadOff * numReadCellPerOperationNeuro);

				// ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                   (0.5 * eReadOff * activityRowRead * numRow * numCol);

				//if (writeBack) {
                            	//	// Calculate Single-Row Read-Out Energy
                            	//	// Parameters: (columnResistance, numColMuxed, numRead)
                            	//	currentSenseAmp.CalculatePower(columnResistance, 1);

                            	//	// Scale sensing energy for the entire subarray (row-by-row)
                            	//	double senseEnergyTotal = currentSenseAmp.readDynamicEnergy * numRow * numCol * bitsPerCell;

                            	//	// 1. Peripheral Energy to drive all lines
                            	//	double peripheralWriteEnergy = wblDecoder.writeDynamicEnergy + wblSwitchMatrix.writeDynamicEnergy +
                            	//	                               wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.writeDynamicEnergy +
                            	//	                               wplDecoder.writeDynamicEnergy + wplSwitchMatrix.writeDynamicEnergy +
                            	//	                               sslDecoder.writeDynamicEnergy + sslSwitchMatrix.writeDynamicEnergy;

				//	// 2. FeCap Flipping Energy
				//	//eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));
				//	double Qsw    = param->twoPr * param->cellAreaFE;
				//	eWritePerCell = Qsw * cell.writeVoltage + 0.5 * CapFe * pow(cell.writeVoltage, 2);

                            	//	// We must rewrite the entire subarray to refresh the states
                            	//	 // double totalCellFlippingEnergy = eWritePerCell * numCol * numRow * param->bitsPerCell;

                            	//	// 3. Total Refresh Energy
                            	//	// (Multiply peripheral energy by numRow because we write row-by-row)
                            	//	// double writeBackEnergy = (peripheralWriteEnergy * numRow * param->bitsPerCell) + totalCellFlippingEnergy;
                            	//	double writeBackEnergy = (peripheralWriteEnergy * numRow * bitsPerCell) ;

                            	//	readDynamicEnergyArray += writeBackEnergy;
				//}


			} else if (cell.mem_rdo == Type::dro) {

				double alphaRd = 1.0;
				double V_RBL   = 0.5;
				double tREAD   = 10e-9;

				dielectric = CapFe * pow(cell.readVoltage, 2);

				e_WBL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_WBL_OFF = 0.0;

				e_RBL_ON  = V_RBL * iRead1 * tREAD;
				e_RBL_OFF = V_RBL * iRead0 * tREAD;

				eReadOn  = e_WBL_ON  + e_RBL_ON  + dielectric;
				eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;

				// // 1. WBL Energy: Driven by the input drivers (DACs/PWM)
                                // // V_read is applied to the WBL of the target cells

                                // double iFlip_ON = 20e-6; // 20 uA instantaneous switching current
                                // double iFlip_OFF = 1.5e-6;
                                // double tFlip = 2e-9;

                                // double V_RBL = 0.5;
                                // double iRBL_ON = 2e-6;
                                // double iRBL_OFF = 10e-9;
                                // double tREAD = 10e-9;

				// dielectric = 0.5 * CapFe * pow(cell.readVoltage,2);

				// e_WBL_ON = cell.readVoltage * iFlip_ON * tFlip;
                                // e_WBL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

                                // // 2. RBL Energy: Driven by a small bias voltage during read
                                // // A small voltage is simultaneously applied to the RBL
                                // e_RBL_ON = V_RBL * iRBL_ON * tREAD;
                                // e_RBL_OFF = V_RBL * iRBL_OFF * tREAD;

                                // eReadOn = e_WBL_ON + e_RBL_ON + dielectric;
                                // eReadOff = e_WBL_OFF + e_RBL_OFF + dielectric;

                                // //double ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow) +
                                // //                         (0.5 * eReadOff * numReadCellPerOperationNeuro);

				// ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                   (0.5 * eReadOff * activityRowRead * numRow * numCol);

				// WriteBack
				// Calculate Single-Row Read-Out Energy
                                // Parameters: (columnResistance, numColMuxed, numRead)
                                currentSenseAmp.CalculatePower(columnResistance, 1);

                                // Scale sensing energy for the entire subarray (row-by-row)
                                double senseEnergyTotal = currentSenseAmp.readDynamicEnergy * numRow * bitsPerCell;

                                // 1. Peripheral Energy to drive all lines
                                double voltageMultiplier = pow(cell.writeVoltage / cell.readVoltage, 2);

				double peripheralWriteEnergy = (wblDecoder.writeDynamicEnergy + wblSwitchMatrix.readDynamicEnergy +
                                                               wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.readDynamicEnergy +
                                                               wplDecoder.writeDynamicEnergy + wplSwitchMatrix.readDynamicEnergy +
                                                               sslDecoder.writeDynamicEnergy + sslSwitchMatrix.readDynamicEnergy) * voltageMultiplier;
				
				// double peripheralWriteEnergy = wblDecoder.writeDynamicEnergy + wblSwitchMatrix.writeDynamicEnergy +
                                //                                wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.writeDynamicEnergy +
                                //                                wplDecoder.writeDynamicEnergy + wplSwitchMatrix.writeDynamicEnergy +
                                //                                sslDecoder.writeDynamicEnergy + sslSwitchMatrix.writeDynamicEnergy;

				// 2. FeCap Flipping Energy
                                // eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));
				double Qsw    = param->twoPr * param->cellAreaFE;
				eWritePerCell = Qsw * cell.writeVoltage + 0.5 * CapFe * pow(cell.writeVoltage, 2);
                                

                                // 3. Total Refresh Energy
                                // double writeBackEnergy = (peripheralWriteEnergy * numRow * param->bitsPerCell) + totalCellFlippingEnergy;
//                                double writeBackEnergy = (peripheralWriteEnergy * numRow * bitsPerCell);

  //                              readDynamicEnergyArray += writeBackEnergy;


			}


                	//readDynamicEnergyArray += ReadEnergyArray;

                	//core read energy 
			readSenseEnergyCore = ReadEnergyArray;

				currentSenseAmp.CalculatePower(columnResistance, 1);
				senseEnergyCore = currentSenseAmp.readDynamicEnergy;
				// readEnergySense = readEnergyCore + senseEnergyCore;

			
			double restoreDuty = 0.0;
			if      (cell.mem_rdo == Type::dro)   restoreDuty = 1.0;
			else if (cell.mem_rdo == Type::qndro) restoreDuty = 1.0 / param->qndroRefreshInterval;

			if (restoreDuty > 0) {
			    double QswR    = param->twoPr * param->cellAreaFE;
			    
			    double eCellR  = QswR * cell.writeVoltage
			                   + CapFe * pow(cell.writeVoltage, 2);
			    double switchR = eCellR * numWriteCellPerOperationNeuro
			                     * numRow * activityRowWrite * param->dataOnesRead;

			    double vInhR   = cell.writeVoltage / 2.0;
			    double capFeHR = CapFe;
			    
			    // double eCellR  = QswR * cell.writeVoltage
			    //                + 0.5 * CapFe * pow(cell.writeVoltage, 2);
			    // double switchR = eCellR * numWriteCellPerOperationNeuro
			    //                  * numRow * activityRowWrite * 0.5;          /* alpha_Fe */

			    // double vInhR   = cell.writeVoltage / 3.0;
			    // double capFeHR = 0.5e-15;
			    int    nUnselR = param->numColSubArray * ((int)param->bitsPerCell - 1)
			                   + param->numColSubArray * (param->numRowSubArrayPhysical - 1);
			    double inhibR  = 0.5 * capFeHR * vInhR * vInhR * nUnselR
			                     * numWriteOperationPerRow * numRow * activityRowWrite;

			    readRestoreEnergyCore = restoreDuty * (switchR + inhibR);
			}
			readEnergyCore = readSenseEnergyCore + readRestoreEnergyCore;

			// readDynamicEnergyArray += readRestoreEnergyCore;   /* REVERTED: same reason.
			//     readDynamicEnergyArray feeds readDynamicEnergy, which ProcessingUnit
			//     consumes.  The core/full tier bookkeeping is done in mem_main.cpp. */

			// double restoreDuty = 0.0;
			// if      (cell.mem_rdo == Type::dro)   restoreDuty = 1.0;
			// else if (cell.mem_rdo == Type::qndro) restoreDuty = 1.0 / param->qndroRefreshInterval;
			// 

			// if (restoreDuty > 0) {
			//     // a restore is a full page write
			//     double eWritePerCellR = (cell.writeVoltage * iFlip_ON_W * tWrite)
			//                           + (0.5 * CapFe * pow(cell.writeVoltage, 2));
			//     double switchR = eWritePerCellR * numWriteCellPerOperationNeuro
			//                      * numRow * activityRowRead * 0.5;            /* alpha_Fe */

			//     double vInhR   = cell.writeVoltage / 3.0;
			//     double capFeHR = 0.5e-15;
			//     int    nTotalR = param->numRowSubArrayPhysical * param->numColSubArray
			//                      * (int)param->bitsPerCell;                   
			//     int    nUnselR = nTotalR - numWriteCellPerOperationNeuro;
			//     double inhibR  = 0.5 * capFeHR * vInhR * vInhR * nUnselR
			//                      * numRow * activityRowRead;

			//     readRestoreEnergyCore = restoreDuty * (switchR + inhibR);
			// }
			// readEnergyCore = readSenseEnergyCore + readRestoreEnergyCore;

			readDynamicEnergyArray += ReadEnergyArray;
			
			//Total Read Energy
                	readDynamicEnergy = 0;
                	readDynamicEnergy += ((numAdd > 1) ? (adder.readDynamicEnergy + dff.readDynamicEnergy) : 0);
                	readDynamicEnergy += rslDecoder.readDynamicEnergy + rslSwitchMatrix.readDynamicEnergy;
                	//readDynamicEnergy += rblDecoder.readDynamicEnergy + rblSwitchMatrix.readDynamicEnergy; // Is Decoder Required here
                	readDynamicEnergy += rblSwitchMatrix.readDynamicEnergy; 
                	readDynamicEnergy += wblDecoder.readDynamicEnergy + wblSwitchMatrix.readDynamicEnergy;
                	readDynamicEnergy += readDynamicEnergyArray;

			
			// if (cell.mem_rdo == Type::ndro) {
			// 	
			// 	// 1. Calculate DC Read Power for a single ON cell
        		// 	// Power = V^2 / R
        		// 	double pReadOn = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOn;

        		// 	// 2. Calculate DC Read Power for a single OFF cell (Leakage during read)
        		// 	double pReadOff = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOff;

        		// 	// 3. Calculate Energy per cell for the duration of the read pulse
        		// 	// Energy = Power * Time
        		// 	double eReadOn = pReadOn * cell.readPulseWidth;
        		// 	double eReadOff = pReadOff * cell.readPulseWidth;

        		// 	// 4. Calculate total number of cells being actively read in this cycle
        		// 	//double numActiveReadCells = numCol * activityColRead;
        		// 	double numActiveReadCells = numReadCellPerOperationNeuro;

        		// 	// 5. Total Array DC Read Energy
        		// 	// Assuming a standard random data distribution (50% ON, 50% OFF)
        		// 	double arrayDCReadEnergy = (0.5 * eReadOn * numActiveReadCells) +
        		// 	                           (0.5 * eReadOff * numActiveReadCells);

        		// 	// 6. Add the cell DC energy to the total array dynamic energy
        		// 	readDynamicEnergyArray += arrayDCReadEnergy;

			//  } else if (cell.mem_rdo == Type::qndro) {

			//  	// 1. Normal DC Read Power
			//  	double pReadOn = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOn;
			//  	double pReadOff = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOff;
			//  	
			//  	double eReadOn = pReadOn * cell.readPulseWidth;
			//  	double eReadOff = pReadOff * cell.readPulseWidth;
			//  	
			//  	// 2. QNDR Domain Flipping Energy
			//  	double vRead_QNDR = cell.readVoltage;
			//  	double iFlip = 20e-6; // 20 uA instantaneous switching current
			//  	// double tFlip = 5e-9;  // 5ns switching time
			//  	
			//  	// Energy to flip twice (Forward + Reverse Self-Restore)
			//  	double eFlipTwice = 2 * (vRead_QNDR * iFlip * cell.readPulseWidth);
			//  	
			//  	double numActiveReadCells = numReadCellPerOperationNeuro;

			//  	double readDisturb = 10;

			//  	double eWriteOne = (cell.writeVoltage / cell.resMemCellOn) * cell.writeVoltage * cell.writePulseWidth;

			//  	double eRefresh = eWriteOne / readDisturb;
			//  	
			//  	// 3. Total Array Energy = DC Read Energy + Double Flipping Energy
			//  	double arrayDCReadEnergy = (0.5 * (eReadOn + eFlipTwice + eRefresh) * numActiveReadCells) +
			//  	                           (0.5 * (eReadOff + eRefresh) * numActiveReadCells);
			//  	
			//  	readDynamicEnergyArray += arrayDCReadEnergy;

			//  } else if (cell.mem_rdo == Type::dro) {
    			//  	
			// 	// Normal DC Read Power
                        //         double pReadOn = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOn;
                        //         double pReadOff = (cell.readVoltage * cell.readVoltage) / cell.resMemCellOff;

                        //         double eReadOn = pReadOn * cell.readPulseWidth;
                        //         double eReadOff = pReadOff * cell.readPulseWidth;

			// 	 
			// 	 // 1. Full Destructive Read Flip Energy
    			//  	double iFlipFull = 20e-6; // Full switching current
    			//  	double tFlipFull = 10e-9; // Full switching time
    			//  	double eFlipDestructive = cell.readVoltage * iFlipFull * tFlipFull;

    			//  	// 2. Write-Back Peripheral Energy (WWL, WPL, WBL, SSL activated)
    			//  	double writeBackPeripheralEnergy = wblDecoder.writeDynamicEnergy + wblSwitchMatrix.writeDynamicEnergy +
    			//  	                                   wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.writeDynamicEnergy +
    			//  	                                   wplDecoder.writeDynamicEnergy + wplSwitchMatrix.writeDynamicEnergy +
    			//  	                                   sslDecoder.writeDynamicEnergy + sslSwitchMatrix.writeDynamicEnergy;

    			//  	// 3. Write-Back Cell Flip Energy
    			//  	double eFlipWriteBack = cell.writeVoltage * iFlipFull * tFlipFull;

    			//  	// 4. Calculate Total Energy
    			//  	double numActiveReadCells = numReadCellPerOperationNeuro;
    			//  	
    			//  	// Assume 50% of cells were '1' and got destroyed, requiring the full read-flip + write-back penalty
    			//  	double arrayDROEnergy = 0.5 * (eReadOn + eFlipDestructive + writeBackPeripheralEnergy + eFlipWriteBack) * numActiveReadCells + 
    			//  	                        0.5 * eReadOff * numActiveReadCells;
    			//  	                        
    			//  	readDynamicEnergyArray += arrayDROEnergy;

			//  }


                	//Peripheral Energy
                	readDynamicEnergy += ((numColMuxed > 1) ? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy) : 0);
                	readDynamicEnergy += (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy));
                	readDynamicEnergy += ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0 )
				+ ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0);


                	//Calculate Write Dynamic Energy

                	//Array Energy: Energy to charge the Internal Capacitors
                	//Energy = 0.5 * C * V^2
                	//From RRAM
                	double totalCapPerCell = cell.numCapacitors * cell.capacitance;
                	//writeDynamicEnergyArray = totalCapPerCell * cell.writeVoltage * cell.writeVoltage;

                	//No. of cells being written
                	//writeDynamicEnergyArray *= numWriteOperationPerRow * numRow * activityRowWrite;

                	//Total Write Energy
                	writeDynamicEnergy = 0;
                	//DRiver Energy (4 Write Matrices)
                	//writeDynamicEnergy += wblDecoder.writeDynamicEnergy + wblSwitchMatrix.writeDynamicEnergy;
                	//writeDynamicEnergy += sslDecoder.writeDynamicEnergy + sslSwitchMatrix.writeDynamicEnergy;
                	//writeDynamicEnergy += wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.writeDynamicEnergy;
                	//writeDynamicEnergy += wplDecoder.writeDynamicEnergy + wplSwitchMatrix.writeDynamicEnergy;

                	double writeErase = 0;
			double alpha_Fe = 0.5;

			double ArrayWriteEnergy = 0;

			//eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));
			double Qsw        = param->twoPr * param->cellAreaFE;          // add twoPr to Param (C/m^2, e.g. 0.20)
			//eWritePerCell     = Qsw * cell.writeVoltage + 0.5 * CapFe * pow(cell.writeVoltage, 2);
			eWritePerCell     = Qsw * cell.writeVoltage + CapFe * pow(cell.writeVoltage, 2);


			//ArrayWriteEnergy = eWritePerCell * numWriteOperationPerRow * numRow * alpha_Fe;
			ArrayWriteEnergy = eWritePerCell * numWriteCellPerOperationNeuro * numRow * activityRowWrite * alpha_Fe;

                	 //writeErase = 0.5 * numRow * numCol * 8 * (cell.writeVoltage * iFlip_ON * tWrite);
			
			
			double peripheralWriteEnergy = wblDecoder.writeDynamicEnergy + wblSwitchMatrix.writeDynamicEnergy +
                                       wwlDecoder.writeDynamicEnergy + wwlSwitchMatrix.writeDynamicEnergy +
                                       wplDecoder.writeDynamicEnergy + wplSwitchMatrix.writeDynamicEnergy +
				       sslDecoder.writeDynamicEnergy + sslSwitchMatrix.writeDynamicEnergy;

        		// // 2. Inputs for FeCap Flipping
        		// double vWrite = cell.writeVoltage;
        		// double tWritePulse = cell.writePulseWidth;
        		// 
        		// // 3. Instantaneous Flipping Current
        		// // 20 microamps per cell during the flip
        		// double iFlip = 20e-6; 

        		// // 4. Manual Energy Calculation for a Single FeCap Flip
        		// // E = V * I * t
        		// double eFlipPerCell = vWrite * iFlip * tWritePulse;

        		// // 5. Total Array Cell Flipping Energy
        		// double totalCellFlippingEnergy = eFlipPerCell * numWriteCellPerOperationNeuro;

			
			// Inhibition scheme (V/3 scheme is typical for cross-point-like sharing)
			//double vInhibit = cell.writeVoltage / 3.0; 
			
			double vInhibit = cell.writeVoltage / 2.0;
			{   
			    double eInhField = (vInhibit / param->tFE) / 1e8;   /* MV/cm */
			    static bool warnedInh2T = false;
			    if (eInhField > param->ecFE && !warnedInh2T) {
			        warnedInh2T = true;
			        cout << "\n INHIBIT DISTURB (2T-nC): half-select field "
			             << eInhField << " MV/cm exceeds Ec = " << param->ecFE
			             << " MV/cm.\n    Lower writeVoltage below "
			             << 2*param->ecFE*param->tFE*1e8 << " V, or thicken tFE.\n\n";
			    }
			}


			int numTotalCapacitors = param->numRowSubArrayPhysical * param->numColSubArray * (int)param->bitsPerCell;
			//int numUnselectedCells = numTotalCapacitors - numWriteCellPerOperationNeuro;
			
			int numUnselectedCells = param->numColSubArray * ((int)param->bitsPerCell - 1)
			                       + param->numColSubArray * (param->numRowSubArrayPhysical - 1);

			//double capFeH = 0.5e-15; 
			double capFeH = CapFe;   
			double inhibitionEnergy = 0.5 * capFeH * (vInhibit * vInhibit) * numUnselectedCells;
			inhibitionEnergy *= (numWriteOperationPerRow * numRow * activityRowWrite);

			// Multiply by the number of write cycles (2 cycles for a parallel write scheme)
			// writeCyclesPerCell = 2;
			// inhibitionEnergy *= writeCyclesPerCell;
			
			// 6. Aggregate Total Write Energy (Adding inhibition penalty)
			// writeDynamicEnergyArray = peripheralWriteEnergy + ArrayWriteEnergy + inhibitionEnergy;
			writeCellEnergyCore  = ArrayWriteEnergy;
			inhibitionEnergyCore = inhibitionEnergy;
			writeEnergyCore      = ArrayWriteEnergy + inhibitionEnergy;
			writeDynamicEnergyArray = peripheralWriteEnergy + writeEnergyCore;
			
			//Array Energy
                	writeDynamicEnergy += writeDynamicEnergyArray;


                	//Calculate Leakage Power
                	leakage = 0;
                	//Leakage from 6 matrices
                	leakage += wblDecoder.leakage + wblSwitchMatrix.leakage;
                	leakage += sslDecoder.leakage + sslSwitchMatrix.leakage;
                	leakage += rslDecoder.leakage + rslSwitchMatrix.leakage;
                	leakage += wwlDecoder.leakage + wwlSwitchMatrix.leakage;
                	leakage += wplDecoder.leakage + wplSwitchMatrix.leakage;
                	leakage += rblDecoder.leakage + rblSwitchMatrix.leakage;

                	//Peripheral Leakage
                	leakage += ((numColMuxed > 1) ? (mux.leakage + muxDecoder.leakage) : 0);
                	leakage += (param->SARADC ? sarADC.leakage : (multilevelSenseAmp.leakage + multilevelSAEncoder.leakage));
                	leakage += ((numCellPerSynapse > 1) ? shiftAddWeight.leakage : 0 )
                                 + ((numReadPulse > 1) ? shiftAddInput.leakage : 0);
                	leakage += ((numAdd > 1) ? (adder.leakage  + dff.leakage) : 0);  

                	//Final Calculations
                	readDynamicEnergyADC = readDynamicEnergyArray
                	                                           + (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy));

                	readDynamicEnergyAccum = ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0 )
                                           + ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0);

                	readDynamicEnergyOther = readDynamicEnergy - readDynamicEnergyADC - readDynamicEnergyAccum;

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

			// cout << "\n================ 2TnC Block Energy & Power ================" << endl;
        		// cout << "--- Read Dynamic Energy (J) ---" << endl;
        		// if (rslDecoder.initialized) cout << "RSL Decoder Read Energy:      " << rslDecoder.readDynamicEnergy << " J" << endl;
        		// if (rslSwitchMatrix.initialized) cout << "RSL Switch Matrix Read Energy:" << rslSwitchMatrix.readDynamicEnergy << " J" << endl;
        		// if (rblDecoder.initialized) cout << "RBL Decoder Read Energy:      " << rblDecoder.readDynamicEnergy << " J" << endl;
        		// if (rblSwitchMatrix.initialized) cout << "RBL Switch Matrix Read Energy:" << rblSwitchMatrix.readDynamicEnergy << " J" << endl;

        		// if (numColMuxed > 1) {
        		//     if (muxDecoder.initialized) cout << "MUX Decoder Read Energy:      " << muxDecoder.readDynamicEnergy << " J" << endl;
        		//     if (mux.initialized) cout << "MUX Read Energy:              " << mux.readDynamicEnergy << " J" << endl;
        		// }
        		// if (currentSenseAmp.initialized) {
        		//     cout << "Current Sense Amp Read Energy:" << currentSenseAmp.readDynamicEnergy << " J" << endl;
        		// }
        		// if (param->SARADC) {
        		//     if (sarADC.initialized) cout << "SAR ADC Read Energy:          " << sarADC.readDynamicEnergy << " J" << endl;
        		// } else {
        		//     if (multilevelSenseAmp.initialized) cout << "Multilevel SA Read Energy:    " << multilevelSenseAmp.readDynamicEnergy << " J" << endl;
        		//     if (multilevelSAEncoder.initialized) cout << "Multilevel SA Enc Read Energy:" << multilevelSAEncoder.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numAdd > 1) {
        		//     if (adder.initialized) cout << "Adder Read Energy:            " << adder.readDynamicEnergy << " J" << endl;
        		//     if (dff.initialized) cout << "DFF Read Energy:              " << dff.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numCellPerSynapse > 1) {
        		//     if (shiftAddWeight.initialized) cout << "ShiftAdd Weight Read Energy:  " << shiftAddWeight.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numReadPulse > 1) {
        		//     if (shiftAddInput.initialized) cout << "ShiftAdd Input Read Energy:   " << shiftAddInput.readDynamicEnergy << " J" << endl;
        		// }

        		// cout << "--- Write Dynamic Energy (J) ---" << endl;
        		// if (wwlDecoder.initialized) cout << "WWL Decoder Write Energy:     " << wwlDecoder.writeDynamicEnergy << " J" << endl;
        		// if (wwlSwitchMatrix.initialized) cout << "WWL Switch Matrix Write E:    " << wwlSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		// if (wblDecoder.initialized) cout << "WBL Decoder Write Energy:     " << wblDecoder.writeDynamicEnergy << " J" << endl;
        		// if (wblSwitchMatrix.initialized) cout << "WBL Switch Matrix Write E:    " << wblSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		// if (wplDecoder.initialized) cout << "WPL Decoder Write Energy:     " << wplDecoder.writeDynamicEnergy << " J" << endl;
        		// if (wplSwitchMatrix.initialized) cout << "WPL Switch Matrix Write E:    " << wplSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		// if (sslDecoder.initialized) cout << "SSL Decoder Write Energy:     " << sslDecoder.writeDynamicEnergy << " J" << endl;
        		// if (sslSwitchMatrix.initialized) cout << "SSL Switch Matrix Write E:    " << sslSwitchMatrix.writeDynamicEnergy << " J" << endl;

        		// if (cell.writeVoltage > 1.5) {
        		//     if (wwlLevelShifter.initialized) cout << "WWL Level Shifter Write E:    " << wwlLevelShifter.writeDynamicEnergy << " J" << endl;
        		//     if (wblLevelShifter.initialized) cout << "WBL Level Shifter Write E:    " << wblLevelShifter.writeDynamicEnergy << " J" << endl;
        		//     if (wplLevelShifter.initialized) cout << "WPL Level Shifter Write E:    " << wplLevelShifter.writeDynamicEnergy << " J" << endl;
        		//     if (sslLevelShifter.initialized) cout << "SSL Level Shifter Write E:    " << sslLevelShifter.writeDynamicEnergy << " J" << endl;
        		// }

        		// cout << "--- Leakage Power (W) ---" << endl;
        		// if (rslDecoder.initialized) cout << "RSL Decoder Leakage:          " << rslDecoder.leakage << " W" << endl;
        		// if (rslSwitchMatrix.initialized) cout << "RSL Switch Matrix Leakage:    " << rslSwitchMatrix.leakage << " W" << endl;
        		// if (rblDecoder.initialized) cout << "RBL Decoder Leakage:          " << rblDecoder.leakage << " W" << endl;
        		// if (rblSwitchMatrix.initialized) cout << "RBL Switch Matrix Leakage:    " << rblSwitchMatrix.leakage << " W" << endl;
        		// if (wwlDecoder.initialized) cout << "WWL Decoder Leakage:          " << wwlDecoder.leakage << " W" << endl;
        		// if (wblDecoder.initialized) cout << "WBL Decoder Leakage:          " << wblDecoder.leakage << " W" << endl;
        		// if (param->SARADC && sarADC.initialized) cout << "SAR ADC Leakage:              " << sarADC.leakage << " W" << endl;

        		// // --- Destructive Read-Out Write-Back Check ---
        		// if (cell.mem_rdo == Type::dro) {
        		//     cout << "---------------- DRO Penalty (Write-Back Energy) ----------------" << endl;
        		//     if (wwlDecoder.initialized) cout << "[Write-Back] WWL Decoder:     " << wwlDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (wwlSwitchMatrix.initialized) cout << "[Write-Back] WWL SwitchMatrix:" << wwlSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     if (wblDecoder.initialized) cout << "[Write-Back] WBL Decoder:     " << wblDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (wblSwitchMatrix.initialized) cout << "[Write-Back] WBL SwitchMatrix:" << wblSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     if (wplDecoder.initialized) cout << "[Write-Back] WPL Decoder:     " << wplDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (wplSwitchMatrix.initialized) cout << "[Write-Back] WPL SwitchMatrix:" << wplSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     cout << "-----------------------------------------------------------------" << endl;
        		// }
        		// cout << "===========================================================" << endl;

		
		} else if (cell.memCellType == Type::_1TnC) {

			double bitsPerCell = param->bitsPerCell;
                        double numReadOperationPerRow = (numCol > numReadCellPerOperationNeuro) ? numCol / numReadCellPerOperationNeuro : 1;
                        double totalCellFlippingEnergy = 0;

                        //Calculate Operation Counts
                        double numReadCells = (int)ceil((double)numCol/numColMuxed);
                        int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);

			 double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite;


                        //Calculate Driver Power (Switching Energy of wires)
			//
			// 1. Calculate how many times the rows fire during a read
			double totalReadActivations = numRow * activityRowRead;
			double decoderReadActivations = numRow * activityRowRead * numColMuxed;
			double switchMatrixReadActivations = numColMuxed;
			
			// 2. BL (Bitline) - Driven during Read and Write
			// blDecoder.CalculatePower(totalReadActivations, totalWriteActivations);
			// blSwitchMatrix.CalculatePower(totalReadActivations, totalWriteActivations, activityRowRead, activityColWrite);

			double blReadActivations = decoderReadActivations * bitsPerCell;
        		double blSwitchMatrixActivations = numColMuxed * bitsPerCell;
			blDecoder.CalculatePower(blReadActivations, totalWriteActivations);
        		blSwitchMatrix.CalculatePower(blSwitchMatrixActivations, totalWriteActivations, activityRowRead, activityColWrite);

			// blDecoder.CalculatePower(decoderReadActivations, totalWriteActivations);
                        // blSwitchMatrix.CalculatePower(switchMatrixReadActivations, totalWriteActivations, activityRowRead, activityColWrite);
			
			// 3. WL (Wordline) - Driven during Read and Write
			// wlDecoder.CalculatePower(totalReadActivations, totalWriteActivations);
			// wlSwitchMatrix.CalculatePower(totalReadActivations, totalWriteActivations, activityRowRead, activityColWrite);

			wlDecoder.CalculatePower(decoderReadActivations, totalWriteActivations);
        		wlSwitchMatrix.CalculatePower(numColMuxed, totalWriteActivations, activityRowRead, activityColWrite);

			// wlDecoder.CalculatePower(decoderReadActivations, totalWriteActivations);
                        // wlSwitchMatrix.CalculatePower(switchMatrixReadActivations, totalWriteActivations, activityRowRead, activityColWrite);

			
			// 4. PL (Plateline) - Passing 1 ensures Leakage is calculated even if it doesn't toggle dynamically during a read
			plDecoder.CalculatePower(1, totalWriteActivations);
			plSwitchMatrix.CalculatePower(1, totalWriteActivations, activityRowRead, activityColWrite);


			
			// double bitsPerCell = param->bitsPerCell;
			// double numReadOperationPerRow = (numCol > numReadCellPerOperationNeuro) ? numCol / numReadCellPerOperationNeuro : 1;
                        // double totalCellFlippingEnergy = 0;

                        // //Calculate Operation Counts
                        // //How many times it loops to process a single logical row
                        // double numReadCells = (int)ceil((double)numCol/numColMuxed);
                        // int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);

			// //Calculate Driver Power (Switching Energy of wires)

                        // // 2-Cycle parallel write requires the peripheral drivers to fire exactly twice per cell write
                        // //int writeCyclesPerCell = 2;
                        // // 1. Calculate Global Operation Counts
			// double totalReadActivations = numRow * activityRowRead;


			// // double numReadCells = (int)ceil((double)numCol/numColMuxed);
                        // // int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);

                        // //Calculate Driver Power (Switching Energy of wires)


                        // // allel write requires the peripheral drivers to fire exactly twice per cell write
                        // //int writeCyclesPerCell = 2;
                        // //    double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite;


                        // //    // PL
                        // //    plDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, numRow*activityRowWrite);
                        // //    plSwitchMatrix.CalculatePower(numColMuxed, 0, activityRowRead, activityColWrite);

                        // //    // BL (Bitline - Row Plane)
                        // //    blDecoder.CalculatePower(0, totalWriteActivations);
                        // //    blSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

			// //    // WL
                        // //    wlDecoder.CalculatePower(0, totalWriteActivations);
                        // //    wlSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

			// // // int writeCyclesPerCell = 2; 
			// double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite; // * writeCyclesPerCell;

			// // DRIVER POWER
			// 
                        //  // plDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, numRow*activityRowWrite);
                        //  // plSwitchMatrix.CalculatePower(numColMuxed, 0, activityRowRead, activityColWrite);

                        //  // blDecoder.CalculatePower(0, totalWriteActivations);
                        //  // blSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

                        //  // wlDecoder.CalculatePower(0, totalWriteActivations);
                        //  // wlSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);

			//   // BL (Plane - Row): Precharged during read, driven during write
			//   blDecoder.CalculatePower(0, totalWriteActivations);
			//   blSwitchMatrix.CalculatePower(0, totalWriteActivations, activityRowRead, activityColWrite);
			//   
			//   // WL (Row): Transistor Gate - Active during BOTH Read and Write
			//   wlDecoder.CalculatePower(totalReadActivations, totalWriteActivations);
			//   wlSwitchMatrix.CalculatePower(totalReadActivations, totalWriteActivations, activityRowRead, activityColWrite);
			//   
			//   // PL (Column): Write Plate Line - Active during Write
			//   plDecoder.CalculatePower(totalReadActivations, totalWriteActivations);
			//   plSwitchMatrix.CalculatePower(totalReadActivations, totalWriteActivations, activityRowRead, activityColWrite);
			
			
			// PERIPHERALS (Mux, ADC, etc)
			if (numColMuxed > 1) {
				mux.CalculatePower(numColMuxed);
				muxDecoder.CalculatePower(numColMuxed, 1);
			}
			
			if (param->SARADC) {
				sarADC.CalculatePower(columnResistance, numAdd);
			} else {
				multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
				multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
			}
			
			
			// ACCUMULATION LOGIC
			if (numCellPerSynapse > 1) {
				shiftAddWeight.CalculatePower(numColMuxed);
			}
			if (numReadPulse > 1) {
				shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));
			}
			if (numAdd > 1) {
				adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed));
				dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
			}
			
                        double resCellWL = resRow;

                        double resCellBL = resRow;

                        double resCellPL = resCol;

                        double resReadTransistor = cell.resMemCellOn; // Assuming worst-case ON resistance

                        double e_WL_ON = 0;
                        double e_WL_OFF = 0;
                        double ReadEnergyArray = 0;
                        double dielectric = 0;

                        double e_BL_ON = 0;
                        double e_BL_OFF = 0;
                        double eReadOn = 0;
                        double eReadOff = 0;

                        double iFlip_ON_W = 20e-6;
                        double tWrite = 10e-9;
                        //double CapFe = 0.5e-15;

			double CapFe = EPS0_VACUUM * param->epsFE * param->cellAreaFE / param->tFE;
			
			double eWritePerCell = 0;
			double writeBackEnergy = 0;
			double senseEnergyTotal = 0;
                        double peripheralWriteEnergy = 0;


                        if (cell.mem_rdo == Type::ndro) {

                                //// 1. WBL Energy: Driven by the input drivers (DACs/PWM)
                                //// V_read is applied to the WBL of the target cells

                                //double iFlip_ON = 2e-6; // 2 uA instantaneous switching current
                                //double iFlip_OFF = 1.5e-6;
                                //double tFlip = 0.7e-9;

                                //dielectric = 0.5 * CapFe * cell.readVoltage * cell.readVoltage;

                                //// 1. BL:
                                //e_BL_ON = cell.readVoltage * iFlip_ON * tFlip;
                                //e_BL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

                                //eReadOn = e_BL_ON + dielectric;
                                //eReadOff = e_BL_OFF + dielectric;

				////ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow) + (0.5 * eReadOff * activityRowRead * numRow);
				// ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                  (0.5 * eReadOff * activityRowRead * numRow * numCol);
				//

				double alphaRd = 0.0;
				double QswRd   = param->twoPr * param->cellAreaFE;
				double p1Rd    = param->dataOnesRead;

                                dielectric = CapFe * cell.readVoltage * cell.readVoltage;

				e_BL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_BL_OFF = 0.0;

                                eReadOn  = e_BL_ON  + dielectric;
                                eReadOff = e_BL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;


                        } else if (cell.mem_rdo == Type::qndro) {

                                // // 1. BL Energy: Driven by the input drivers (DACs/PWM)
                                // // V_read is applied to the WBL of the target cells

                                // double iFlip_ON = 10e-6; // 20 uA instantaneous switching current
                                // double iFlip_OFF = 1.5e-6;
                                // double tFlip = 2e-9;

                                // dielectric = 0.5 * CapFe * cell.readVoltage * cell.readVoltage;

                                // // 1. BL Energy:
                                // e_BL_ON = cell.readVoltage * iFlip_ON * tFlip;
                                // e_BL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

                                // eReadOn = e_BL_ON + dielectric;
                                // eReadOff = e_BL_OFF + dielectric;

                                // // ReadEnergyArray = (0.5 * eReadOn * numReadCellPerOperationNeuro) +
                                // //                          (0.5 * eReadOff * numReadCellPerOperationNeuro);

				//  ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                   (0.5 * eReadOff * activityRowRead * numRow * numCol);

				double alphaRd = 0.20;   
				double QswRd   = param->twoPr * param->cellAreaFE;
				double p1Rd    = param->dataOnesRead;

                                dielectric = CapFe * cell.readVoltage * cell.readVoltage;

				e_BL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_BL_OFF = 0.0;

                                eReadOn  = e_BL_ON  + dielectric;
                                eReadOff = e_BL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;

                                //if (writeBack) {
                                //        // Calculate Single-Row Read-Out Energy
                                //        currentSenseAmp.CalculatePower(columnResistance, 1);

                                //        senseEnergyTotal = currentSenseAmp.readDynamicEnergy * numRow * numCol * bitsPerCell;

                                //        // 1. Peripheral Energy to drive all lines
                                //        peripheralWriteEnergy = blDecoder.writeDynamicEnergy + blSwitchMatrix.writeDynamicEnergy +
                                //                                       wlDecoder.writeDynamicEnergy + wlSwitchMatrix.writeDynamicEnergy +
                                //                                       plDecoder.writeDynamicEnergy + plSwitchMatrix.writeDynamicEnergy;

                                //        // 2. FeCap Flipping Energy
                                //        eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));

                                //         // double totalCellFlippingEnergy = eWritePerCell * numCol * numRow * param->bitsPerCell;

                                //        // 3. Total Refresh Energy
                                //        // double writeBackEnergy = (peripheralWriteEnergy * numRow * param->bitsPerCell) + totalCellFlippingEnergy;
                                //        writeBackEnergy = (peripheralWriteEnergy * numRow * bitsPerCell) ;

                                //        readDynamicEnergyArray += writeBackEnergy;
                                //}


                        } else if (cell.mem_rdo == Type::dro) {

                                // // 1. WBL Energy: Driven by the input drivers (DACs/PWM)

                                // double iFlip_ON = 20e-6; // 20 uA instantaneous switching current
                                // double iFlip_OFF = 1.5e-6;
                                // double tFlip = 2e-9;

                                // dielectric = 0.5 * CapFe * pow(cell.readVoltage,2);

                                // e_BL_ON = cell.readVoltage * iFlip_ON * tFlip;
                                // e_BL_OFF = cell.readVoltage * iFlip_OFF * tFlip;

				// eReadOn = e_BL_ON + dielectric;
                                // eReadOff = e_BL_OFF + dielectric;

                                // ReadEnergyArray = (0.5 * eReadOn * activityRowRead * numRow * numCol) +
                                //                   (0.5 * eReadOff * activityRowRead * numRow * numCol);

                                				double alphaRd = 1.0;    /* full reversal: V_read = V_write */
				double QswRd   = param->twoPr * param->cellAreaFE;
				double p1Rd    = param->dataOnesRead;

                                dielectric = CapFe * pow(cell.readVoltage, 2);

				e_BL_ON  = alphaRd * QswRd * cell.readVoltage;
				e_BL_OFF = 0.0;

				eReadOn  = e_BL_ON  + dielectric;
                                eReadOff = e_BL_OFF + dielectric;

				ReadEnergyArray = ((p1Rd * eReadOn) + ((1.0 - p1Rd) * eReadOff))
				                  * activityRowRead * numRow * numCol;
				
				if (activityRowRead > 0) {
				// WriteBack
        			// 1. Sense Energy: currentSenseAmp was initialized for numCol, so 1 activation senses the whole row.
        			// currentSenseAmp.CalculatePower(columnResistance, 1);

        			double totalWriteBackOps = numRow * activityRowRead;

        			// double senseEnergyTotal = currentSenseAmp.readDynamicEnergy * totalWriteBackOps;

        			// 2. Peripheral Energy: Extract the energy of a SINGLE activation, then scale it.
        			double totalReadOps = numRow * activityRowRead * numColMuxed; 

        			double E_decoder_per_op = (blDecoder.readDynamicEnergy + wlDecoder.readDynamicEnergy) / totalReadOps;
        			double E_pl_decoder_per_op = plDecoder.readDynamicEnergy / 1.0; 

        			double E_sw_per_op = (blSwitchMatrix.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy) / numColMuxed;
        			double E_pl_sw_per_op = plSwitchMatrix.readDynamicEnergy / 1.0;

        			double singleReadPeripheralEnergy = E_decoder_per_op + E_pl_decoder_per_op + E_sw_per_op + E_pl_sw_per_op;

        			double voltageMultiplier = pow(cell.writeVoltage / cell.readVoltage, 2);
        			double singleWritePeripheralEnergy = singleReadPeripheralEnergy * voltageMultiplier;

        			peripheralWriteEnergy = singleWritePeripheralEnergy * totalWriteBackOps;

        			// 3. FeCap Flipping Energy
        			eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage, 2));
        			
        			// Total cells written back = total active rows * total columns. Assume 50% need to be flipped back to '1'.
        			double totalCellFlippingEnergy = eWritePerCell * (numRow * activityRowRead * numCol) * 0.5;

        			// 4. Aggregate Write-Back Energy
        			// double writeBackEnergy = peripheralWriteEnergy + senseEnergyTotal + totalCellFlippingEnergy;
        			//writeBackEnergy = peripheralWriteEnergy + totalCellFlippingEnergy;
				}

        			//readDynamicEnergyArray += writeBackEnergy;
				
				// // WriteBack
                                // // Calculate Single-Row Read-Out Energy
                                // currentSenseAmp.CalculatePower(columnResistance, 1);

                                // // senseEnergyTotal = currentSenseAmp.readDynamicEnergy * numRow * bitsPerCell;
                                // senseEnergyTotal = currentSenseAmp.readDynamicEnergy * totalWriteBackOps;

                                // // 1. Peripheral Energy to drive all lines
                                // // double peripheralWriteEnergy = blDecoder.writeDynamicEnergy + blSwitchMatrix.writeDynamicEnergy +
                                // //                                wlDecoder.writeDynamicEnergy + wlSwitchMatrix.writeDynamicEnergy +
                                // //                                plDecoder.writeDy:se nu
				// //                                namicEnergy + plSwitchMatrix.writeDynamicEnergy;

                                // double singleReadPeripheralEnergy = (blDecoder.readDynamicEnergy + blSwitchMatrix.readDynamicEnergy +
        			//                                      wlDecoder.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy) / totalReadOps;
        			// singleReadPeripheralEnergy += (plDecoder.readDynamicEnergy + plSwitchMatrix.readDynamicEnergy) / 1.0; 

        			// double voltageMultiplier = pow(cell.writeVoltage / cell.readVoltage, 2);
        			// peripheralWriteEnergy = singleReadPeripheralEnergy * voltageMultiplier * totalWriteBackOps;

				// // peripheralWriteEnergy = (blDecoder.readDynamicEnergy + blSwitchMatrix.readDynamicEnergy +
				// //                          wlDecoder.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy +
				// //                          plDecoder.readDynamicEnergy + plSwitchMatrix.readDynamicEnergy) * voltageMultiplier;

				// // 2. FeCap Flipping Energy
                                // eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));

				// double totalCellFlippingEnergy = eWritePerCell * totalWriteBackOps * (numCol / numColMuxed) * 0.5;

				// double writeBackEnergy = peripheralWriteEnergy + senseEnergyTotal + totalCellFlippingEnergy;

                                //  // double totalCellFlippingEnergy = esWritePerCell * numCol * numRow * param->bitsPerCell;

                                // // 3. Total Refresh Energy
                                // // double writeBackEnergy = (peripheralWriteEnergy * numRow * param->bitsPerCell) + totalCellFlippingEnergy;
                                // // writeBackEnergy = (peripheralWriteEnergy * numRow * bitsPerCell) + senseEnergyTotal + (eWritePerCell * numRow * bitsPerCell);

                                // readDynamicEnergyArray += writeBackEnergy;


                        }


			readSenseEnergyCore = ReadEnergyArray;

			currentSenseAmp.CalculatePower(columnResistance, 1);
			senseEnergyCore = currentSenseAmp.readDynamicEnergy;

			double restoreDuty = 0.0;
			if      (cell.mem_rdo == Type::dro)   restoreDuty = 1.0;
			else if (cell.mem_rdo == Type::qndro) restoreDuty = 1.0 / param->qndroRefreshInterval;

			if (restoreDuty > 0) {
			    double QswR    = param->twoPr * param->cellAreaFE;
			    //double eCellR  = QswR * cell.writeVoltage
			    //               + 0.5 * CapFe * pow(cell.writeVoltage, 2);
			    //double switchR = eCellR * numWriteCellPerOperationNeuro
			    //                 * numRow * activityRowWrite * 0.5;          /* alpha_Fe */

			    //double vInhR   = cell.writeVoltage / 3.0;
			    //double capFeHR = 0.5e-15;
			    
			    double eCellR  = QswR * cell.writeVoltage
			                   + CapFe * pow(cell.writeVoltage, 2);
			    double switchR = eCellR * numWriteCellPerOperationNeuro
			                     * numRow * activityRowWrite * param->dataOnesRead;

			    double vInhR   = cell.writeVoltage / 2.0;   /* V/2, not V/3 [1] */
			    double capFeHR = CapFe;
			    
			    int    nUnselR = param->numColSubArray * ((int)param->bitsPerCell - 1)
			                   + param->numColSubArray * (param->numRowSubArrayPhysical - 1);
			    double inhibR  = 0.5 * capFeHR * vInhR * vInhR * nUnselR
			                     * numWriteOperationPerRow * numRow * activityRowWrite;

			    readRestoreEnergyCore = restoreDuty * (switchR + inhibR);
			}
			readEnergyCore = readSenseEnergyCore + readRestoreEnergyCore;

                        readDynamicEnergyArray += ReadEnergyArray;

                        //Total Read Energy
                        readDynamicEnergy = 0;
                        readDynamicEnergy += ((numAdd > 1) ? (adder.readDynamicEnergy + dff.readDynamicEnergy) : 0);
                        readDynamicEnergy += wlDecoder.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy;
                        readDynamicEnergy += blDecoder.readDynamicEnergy + blSwitchMatrix.readDynamicEnergy;
                        readDynamicEnergy += readDynamicEnergyArray;


                        //Peripheral Energy
                        readDynamicEnergy += ((numColMuxed > 1) ? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy) : 0);
                        readDynamicEnergy += (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy));
                        readDynamicEnergy += ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0 )
                                + ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0);


                        //Calculate Write Dynamic Energy

                        //Array Energy: Energy to charge the Internal Capacitors
                        //Energy = 0.5 * C * V^2
                        //From RRAM
                        double totalCapPerCell = cell.numCapacitors * cell.capacitance;
                        //writeDynamicEnergyArray = totalCapPerCell * cell.writeVoltage * cell.writeVoltage;

                        //No. of cells being written
                        //writeDynamicEnergyArray *= numWriteOperationPerRow * numRow * activityRowWrite;

                        //Total Write Energy
                        writeDynamicEnergy = 0;

                        double writeErase = 0;
                        double alpha_Fe = 0.5;

                        double ArrayWriteEnergy = 0;

                        eWritePerCell = (cell.writeVoltage * iFlip_ON_W * tWrite) + (0.5 * CapFe * pow(cell.writeVoltage,2));
                        
			double Qsw    = param->twoPr * param->cellAreaFE;
			eWritePerCell = Qsw * cell.writeVoltage + 0.5 * CapFe * pow(cell.writeVoltage, 2);
			
			//ArrayWriteEnergy = eWritePerCell * numWriteOperationPerRow * numRow * alpha_Fe;

                        ArrayWriteEnergy = eWritePerCell * numWriteCellPerOperationNeuro * numRow * activityRowWrite * alpha_Fe;


			//writeErase = 0.5 * numRow * numCol * 8 * (cell.writeVoltage * iFlip_ON * tWrite);

                        peripheralWriteEnergy = blDecoder.writeDynamicEnergy + blSwitchMatrix.writeDynamicEnergy +
                                       wlDecoder.writeDynamicEnergy + wlSwitchMatrix.writeDynamicEnergy +
                                       plDecoder.writeDynamicEnergy + plSwitchMatrix.writeDynamicEnergy;

                        // // 2. Inputs for FeCap Flipping
                        // double vWrite = cell.writeVoltage;
                        // double tWritePulse = cell.writePulseWidth;
                        //
                        // // 3. Instantaneous Flipping Current
                        // // 20 microamps per cell during the flip
                        // double iFlip = 20e-6;

                        // // 4. Manual Energy Calculation for a Single FeCap Flip
                        // // E = V * I * t
                        // double eFlipPerCell = vWrite * iFlip * tWritePulse;

                        // // 5. Total Array Cell Flipping Energy
                        // double totalCellFlippingEnergy = eFlipPerCell * numWriteCellPerOperationNeuro;


                        //double vInhibit = cell.writeVoltage / 3.0;

                        double vInhibit = cell.writeVoltage / 2.0;
			
                        int numTotalCapacitors = param->numRowSubArrayPhysical * param->numColSubArray * (int)param->bitsPerCell;
                        //int numUnselectedCells = numTotalCapacitors - numWriteCellPerOperationNeuro;

			int numUnselectedCells = param->numColSubArray * ((int)param->bitsPerCell - 1)
			                       + param->numColSubArray * (param->numRowSubArrayPhysical - 1);

                        //double capFeH = 0.5e-15;
                                                double capFeH = CapFe;   /* derived, not 0.5e-15 */
			double inhibitionEnergy = 0.5 * capFeH * (vInhibit * vInhibit) * numUnselectedCells;
			inhibitionEnergy *= (numWriteOperationPerRow * numRow * activityRowWrite);

                        // writeCyclesPerCell = 2;
                        // inhibitionEnergy *= writeCyclesPerCell;

                        // 6. Aggregate Total Write Energy (Adding inhibition penalty)
                        //writeDynamicEnergyArray = peripheralWriteEnergy + ArrayWriteEnergy + inhibitionEnergy;
			writeCellEnergyCore  = ArrayWriteEnergy;
			inhibitionEnergyCore = inhibitionEnergy;
			writeEnergyCore      = ArrayWriteEnergy + inhibitionEnergy;
			writeDynamicEnergyArray = peripheralWriteEnergy + writeEnergyCore;

                        //Array Energy
                        writeDynamicEnergy += writeDynamicEnergyArray;


                        //Calculate Leakage Power
                        leakage = 0;
                        //Leakage from 3 matrices
                        leakage += blDecoder.leakage + blSwitchMatrix.leakage;
                        leakage += wlDecoder.leakage + wlSwitchMatrix.leakage;
                        leakage += plDecoder.leakage + plSwitchMatrix.leakage;

                        //Peripheral Leakage
                        leakage += ((numColMuxed > 1) ? (mux.leakage + muxDecoder.leakage) : 0);
                        leakage += (param->SARADC ? sarADC.leakage : (multilevelSenseAmp.leakage + multilevelSAEncoder.leakage));
                        leakage += ((numCellPerSynapse > 1) ? shiftAddWeight.leakage : 0 )
                                 + ((numReadPulse > 1) ? shiftAddInput.leakage : 0);
                        leakage += ((numAdd > 1) ? (adder.leakage  + dff.leakage) : 0);

                        //Final Calculations
                        readDynamicEnergyADC = readDynamicEnergyArray
                                                                   + (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy));

                        readDynamicEnergyAccum = ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0 )
                                           + ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0);

                        readDynamicEnergyOther = readDynamicEnergy - readDynamicEnergyADC - readDynamicEnergyAccum;

			// cout << "\n================ 1TnC Block Energy & Power ================" << endl;
        		// cout << "--- Read Dynamic Energy (J) ---" << endl;
        		// if (wlDecoder.initialized) cout << "WL Decoder Read Energy:       " << wlDecoder.readDynamicEnergy << " J" << endl;
        		// if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Read Energy: " << wlSwitchMatrix.readDynamicEnergy << " J" << endl;
        		// if (blDecoder.initialized) cout << "BL Decoder Read Energy:       " << blDecoder.readDynamicEnergy << " J" << endl;
        		// if (blSwitchMatrix.initialized) cout << "BL Switch Matrix Read Energy: " << blSwitchMatrix.readDynamicEnergy << " J" << endl;
        		// if (plDecoder.initialized) cout << "PL Decoder Read Energy:       " << plDecoder.readDynamicEnergy << " J" << endl;
        		// if (plSwitchMatrix.initialized) cout << "PL Switch Matrix Read Energy: " << plSwitchMatrix.readDynamicEnergy << " J" << endl;

        		// if (numColMuxed > 1) {
        		//     if (muxDecoder.initialized) cout << "MUX Decoder Read Energy:      " << muxDecoder.readDynamicEnergy << " J" << endl;
        		//     if (mux.initialized) cout << "MUX Read Energy:              " << mux.readDynamicEnergy << " J" << endl;
        		// }
        		// if (currentSenseAmp.initialized) {
        		//     cout << "Current Sense Amp Read Energy:" << currentSenseAmp.readDynamicEnergy << " J" << endl;
        		// }
        		// if (param->SARADC) {
        		//     if (sarADC.initialized) cout << "SAR ADC Read Energy:          " << sarADC.readDynamicEnergy << " J" << endl;
        		// } else {
        		//     if (multilevelSenseAmp.initialized) cout << "Multilevel SA Read Energy:    " << multilevelSenseAmp.readDynamicEnergy << " J" << endl;
        		//     if (multilevelSAEncoder.initialized) cout << "Multilevel SA Enc Read Energy:" << multilevelSAEncoder.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numAdd > 1) {
        		//     if (adder.initialized) cout << "Adder Read Energy:            " << adder.readDynamicEnergy << " J" << endl;
        		//     if (dff.initialized) cout << "DFF Read Energy:              " << dff.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numCellPerSynapse > 1) {
        		//     if (shiftAddWeight.initialized) cout << "ShiftAdd Weight Read Energy:  " << shiftAddWeight.readDynamicEnergy << " J" << endl;
        		// }
        		// if (numReadPulse > 1) {
        		//     if (shiftAddInput.initialized) cout << "ShiftAdd Input Read Energy:   " << shiftAddInput.readDynamicEnergy << " J" << endl;
        		// }

        		// cout << "--- Write Dynamic Energy (J) ---" << endl;
        		// if (wlDecoder.initialized) cout << "WL Decoder Write Energy:      " << wlDecoder.writeDynamicEnergy << " J" << endl;
        		// if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Write Energy:" << wlSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		// if (blDecoder.initialized) cout << "BL Decoder Write Energy:      " << blDecoder.writeDynamicEnergy << " J" << endl;
        		// if (blSwitchMatrix.initialized) cout << "BL Switch Matrix Write Energy:" << blSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		// if (plDecoder.initialized) cout << "PL Decoder Write Energy:      " << plDecoder.writeDynamicEnergy << " J" << endl;
        		// if (plSwitchMatrix.initialized) cout << "PL Switch Matrix Write Energy:" << plSwitchMatrix.writeDynamicEnergy << " J" << endl;

        		// if (cell.writeVoltage > 1.5) {
        		//     if (wlLevelShifter.initialized) cout << "WL Level Shifter Write Energy:" << wlLevelShifter.writeDynamicEnergy << " J" << endl;
        		//     if (blLevelShifter.initialized) cout << "BL Level Shifter Write Energy:" << blLevelShifter.writeDynamicEnergy << " J" << endl;
        		//     if (plLevelShifter.initialized) cout << "PL Level Shifter Write Energy:" << plLevelShifter.writeDynamicEnergy << " J" << endl;
        		// }

        		// cout << "--- Leakage Power (W) ---" << endl;
        		// if (wlDecoder.initialized) cout << "WL Decoder Leakage:           " << wlDecoder.leakage << " W" << endl;
        		// if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Leakage:     " << wlSwitchMatrix.leakage << " W" << endl;
        		// if (blDecoder.initialized) cout << "BL Decoder Leakage:           " << blDecoder.leakage << " W" << endl;
        		// if (blSwitchMatrix.initialized) cout << "BL Switch Matrix Leakage:     " << blSwitchMatrix.leakage << " W" << endl;
        		// if (plDecoder.initialized) cout << "PL Decoder Leakage:           " << plDecoder.leakage << " W" << endl;
        		// if (plSwitchMatrix.initialized) cout << "PL Switch Matrix Leakage:     " << plSwitchMatrix.leakage << " W" << endl;
        		// if (numColMuxed > 1) {
        		//     if (muxDecoder.initialized) cout << "MUX Decoder Leakage:          " << muxDecoder.leakage << " W" << endl;
        		//     if (mux.initialized) cout << "MUX Leakage:                  " << mux.leakage << " W" << endl;
        		// }
        		// if (currentSenseAmp.initialized) cout << "Current Sense Amp Leakage:    " << currentSenseAmp.leakage << " W" << endl;
        		// if (param->SARADC && sarADC.initialized) cout << "SAR ADC Leakage:              " << sarADC.leakage << " W" << endl;

        		// // --- Destructive Read-Out Write-Back Check ---
        		// if (cell.mem_rdo == Type::dro) {
        		//     cout << "---------------- DRO Penalty (Write-Back Energy) ----------------" << endl;
        		//     if (wlDecoder.initialized) cout << "[Write-Back] WL Decoder:      " << wlDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (wlSwitchMatrix.initialized) cout << "[Write-Back] WL SwitchMatrix: " << wlSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     if (blDecoder.initialized) cout << "[Write-Back] BL Decoder:      " << blDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (blSwitchMatrix.initialized) cout << "[Write-Back] BL SwitchMatrix: " << blSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     if (plDecoder.initialized) cout << "[Write-Back] PL Decoder:      " << plDecoder.writeDynamicEnergy << " J" << endl;
        		//     if (plSwitchMatrix.initialized) cout << "[Write-Back] PL SwitchMatrix: " << plSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		//     cout << "-----------------------------------------------------------------" << endl;
        		// }
        		// cout << "===========================================================" << endl;

		
		} else if (cell.memCellType == Type::_1T1C) {
        		int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
        		double totalWriteActivations = numWriteOperationPerRow * numRow * activityRowWrite;
        		double totalReadActivations = numRow * activityRowRead;
        		double decoderReadActivations = numRow * activityRowRead * numColMuxed;

        		// BL (Col Plane)
        		blDecoder.CalculatePower(decoderReadActivations, totalWriteActivations);
        		blSwitchMatrix.CalculatePower(numColMuxed, totalWriteActivations, activityRowRead, activityColWrite);

        		// WL (Row Plane)
        		wlDecoder.CalculatePower(decoderReadActivations, totalWriteActivations);
        		wlSwitchMatrix.CalculatePower(numColMuxed, totalWriteActivations, activityRowRead, activityColWrite);

        		if (numColMuxed > 1) {
        		    mux.CalculatePower(numColMuxed);
        		    muxDecoder.CalculatePower(numColMuxed, 1);
        		}

        		if (param->SARADC) sarADC.CalculatePower(columnResistance, numAdd);
        		else {
        		    multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
        		    multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
        		}

        		if (numCellPerSynapse > 1) shiftAddWeight.CalculatePower(numColMuxed);
        		if (numReadPulse > 1) shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));
        		if (numAdd > 1) {
        		    adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed));
        		    dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
        		}

        		
			// DRAM Destructive Read Write-Back Energy
        		double totalWriteBackOps = numRow * activityRowRead;
        		
        		// Prevent division by zero if activity is 0
        		double E_decoder_per_op = 0;
        		if (totalWriteBackOps > 0) {
        		    E_decoder_per_op = (blDecoder.readDynamicEnergy + wlDecoder.readDynamicEnergy) / (totalWriteBackOps * numColMuxed);
        		}
        		double E_sw_per_op = (blSwitchMatrix.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy) / numColMuxed;
        		
        		// Scale read routing energy to write voltages
        		double voltageMultiplier = pow(cell.writeVoltage / tech.vdd, 2);
			double peripheralWriteBackEnergy = (E_decoder_per_op + E_sw_per_op) * voltageMultiplier * totalWriteBackOps;

        		double capCell = 5e-15; // 5fF cell capacitance
        		double eWritePerCell = (0.5 * capCell * pow(tech.vdd, 2));
        		
			// Write-back the destroyed array data (assuming 50% 1s)
        		double arrayWriteBackEnergy = eWritePerCell * (numRow * activityRowRead * numCol) * 0.5;

        		double writeBackEnergy = peripheralWriteBackEnergy + arrayWriteBackEnergy;
        		readDynamicEnergyArray = writeBackEnergy;

        		// Calculate total Dynamic Energy
        		readDynamicEnergy = readDynamicEnergyArray
        		                  + wlDecoder.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy
        		                  + blDecoder.readDynamicEnergy + blSwitchMatrix.readDynamicEnergy
        		                  + ((numColMuxed > 1) ? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy) : 0)
        		                  + (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy))
        		                  + ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0) + ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0)
        		                  + ((numAdd > 1) ? (adder.readDynamicEnergy + dff.readDynamicEnergy) : 0);

        		// 3. Standard Write Dynamic Energy
        		// Selected Cell Capacitors: Physically charging the 5fF cells
        		double cellEnergyWrite = (numRow * activityRowWrite * numCol * activityColWrite) * 0.5 * 5e-15 * pow(cell.writeVoltage, 2);


        		writeDynamicEnergyArray = cellEnergyWrite;

        		writeDynamicEnergy = writeDynamicEnergyArray
        		                   + wlDecoder.writeDynamicEnergy + wlSwitchMatrix.writeDynamicEnergy
        		                   + blDecoder.writeDynamicEnergy + blSwitchMatrix.writeDynamicEnergy
        		                   + ((numColMuxed > 1) ? (mux.writeDynamicEnergy + muxDecoder.writeDynamicEnergy) : 0);



        		// double ArrayWriteEnergy = eWritePerCell * numWriteCellPerOperationNeuro * numRow * activityRowWrite * 0.5;
        		// double peripheralWriteEnergy = blDecoder.writeDynamicEnergy + blSwitchMatrix.writeDynamicEnergy + wlDecoder.writeDynamicEnergy + wlSwitchMatrix.writeDynamicEnergy;
        		// double inhibitionEnergy = 0.5 * 0.5e-15 * pow(cell.writeVoltage / 3.0, 2) * ((param->numRowSubArrayPhysical * param->numColSubArray) - numWriteCellPerOperationNeuro) * (numWriteOperationPerRow * numRow * activityRowWrite);

        		// writeDynamicEnergyArray = peripheralWriteEnergy + ArrayWriteEnergy + inhibitionEnergy;
        		// writeDynamicEnergy = writeDynamicEnergyArray;

        		// 4. Leakage Power
        		leakage = blDecoder.leakage + blSwitchMatrix.leakage
        		        + wlDecoder.leakage + wlSwitchMatrix.leakage
        		        + ((numColMuxed > 1) ? (mux.leakage + muxDecoder.leakage) : 0)
        		        + (param->SARADC ? sarADC.leakage : (multilevelSenseAmp.leakage + multilevelSAEncoder.leakage))
        		        + ((numCellPerSynapse > 1) ? shiftAddWeight.leakage : 0) + ((numReadPulse > 1) ? shiftAddInput.leakage : 0)
        		        + ((numAdd > 1) ? (adder.leakage + dff.leakage) : 0);

        		// 5. DRAM Refresh Power Penalty
        		
        		// Total time before a cell loses its charge
        		// double tRefresh = 32e-3; // 32 ms retention time
			double tRefresh = 20.48e-3; // 20.48 ms

        		double tRefreshPerRow = tRefresh / numRow; 

        		double E_decoder_per_row = 0;
        		if (totalWriteBackOps > 0) {
        		    E_decoder_per_row = (blDecoder.readDynamicEnergy + wlDecoder.readDynamicEnergy) / (totalWriteBackOps * numColMuxed);
        		}
        		double E_sw_per_row = (blSwitchMatrix.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy) / numColMuxed;
        		
        		// 1. Energy to Read the row
        		double peripheralReadRowEnergy = E_decoder_per_row + E_sw_per_row;
        		
        		// 2. Energy to Write-Back the row (Scaled to write voltage)
        		double peripheralWriteRowEnergy = peripheralReadRowEnergy * pow(cell.writeVoltage / tech.vdd, 2);
        		
        		// 3. Energy to physically charge the cell capacitors in the row (assuming 50% 1s)
        		double arrayWriteRowEnergy = eWritePerCell * numCol * 0.5; 

        		// Total energy consumed during a single row's refresh cycle
        		double totalRefreshEnergyPerRow = peripheralReadRowEnergy + peripheralWriteRowEnergy + arrayWriteRowEnergy;

        		// Power = Energy / Time 
        		double refreshPower = totalRefreshEnergyPerRow / tRefreshPerRow;
        		
        		leakage += refreshPower;

			readDynamicEnergyADC = (param->SARADC ? sarADC.readDynamicEnergy : (multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy));
        		readDynamicEnergyAccum = ((numCellPerSynapse > 1) ? shiftAddWeight.readDynamicEnergy : 0) + ((numReadPulse > 1) ? shiftAddInput.readDynamicEnergy : 0);
        		readDynamicEnergyOther = readDynamicEnergy - readDynamicEnergyADC - readDynamicEnergyAccum;

        		//     std::cout << "\n================ 1T1C DRAM ENERGY DEBUG =================\n";
        		//     std::cout << "--- 1. Array Dimensions & Physics ---\n";
        		//     std::cout << "numRow: " << numRow << "\n";
        		//     std::cout << "numCol: " << numCol << "\n";
        		//     std::cout << "activityRowRead: " << activityRowRead << "\n";
        		//     std::cout << "tech.vdd: " << tech.vdd << " V\n";
        		//     std::cout << "cell.readVoltage: " << cell.readVoltage << " V\n";
        		//     std::cout << "cell.writeVoltage: " << cell.writeVoltage << " V\n";
        		//     std::cout << "capWL (Total Wordline Cap): " << capWL * 1e15 << " fF\n";
        		//     std::cout << "capBL (Total Bitline Cap): " << capBL * 1e15 << " fF\n";

        		//     //std::cout << "\n--- 2. READ PHASE ENERGY (Per Operation) ---\n";
        		//     //std::cout << "wlEnergyRead: " << wlEnergyRead * 1e12 << " pJ\n";
        		//     //std::cout << "blEnergyRead: " << blEnergyRead * 1e12 << " pJ\n";
        		//     //std::cout << "cellEnergyRead: " << cellEnergyRead * 1e12 << " pJ\n";
        		//     //std::cout << ">> ReadEnergyArrayBase: " << ReadEnergyArrayBase * 1e12 << " pJ\n";

        		//     std::cout << "\n--- 3. WRITE-BACK (DRO) PHASE ENERGY ---\n";
        		//     //std::cout << "wlEnergyWB: " << wlEnergyWB * 1e12 << " pJ\n";
        		//     //std::cout << "blEnergyWB: " << blEnergyWB * 1e12 << " pJ\n";
        		//     //std::cout << "cellEnergyWB: " << cellEnergyWB * 1e12 << " pJ\n";
        		//     std::cout << "eWritePerCell: " << eWritePerCell * 1e12 << " pJ\n";
        		//     std::cout << "arrayWriteBackEnergy (Internal): " << arrayWriteBackEnergy * 1e12 << " pJ\n";
        		//     std::cout << "peripheralWriteBackEnergy (Drivers): " << peripheralWriteBackEnergy * 1e12 << " pJ\n";
        		//     std::cout << ">> writeBackEnergy: " << writeBackEnergy * 1e12 << " pJ\n";

        		//     std::cout << "\n--- 4. AGGREGATED TOTALS ---\n";
        		//     std::cout << "readDynamicEnergyArray (Read + DRO): " << readDynamicEnergyArray * 1e12 << " pJ\n";
        		//     std::cout << "readDynamicEnergyADC (Sense Amps): " << readDynamicEnergyADC * 1e12 << " pJ\n";

        		//     // Prove the math: E = 0.5 * C * V^2
        		//     double manual_bl_energy = 0.5 * capBL * pow(cell.writeVoltage, 2) * numCol;
        		//     std::cout << "\n--- 5. MANUAL MATH CHECK ---\n";
        		//     std::cout << "0.5 * capBL * (V_write)^2 * numCol = " << manual_bl_energy * 1e12 << " pJ\n";
        		//     std::cout << "=========================================================\n\n";

		} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			// Anni update
			leakageSRAMInUse = 0;
			if (conventionalSequential) {
				double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
				double numWriteCells = (int)ceil((double)numCol/*numWriteColMuxed*/); 
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double capBL = lengthCol * 0.2e-15/1e-6;
				wlDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				if (cell.accessType == CMOS_access) {
					wlNewDecoderDriver.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				} else {
					wlDecoderDriver.CalculatePower(numReadCells, numWriteCells, numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				}
				slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}

				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numRow*activityRowRead);
				} else {
					multilevelSenseAmp.CalculatePower(columnResistance, numRow*activityRowRead);
					if (avgWeightBit > 1) {
						multilevelSAEncoder.CalculatePower(numRow*activityRowRead*numColMuxed);
					}
				}
				adder.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells);
				// Anni update
				dff.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells*(ceil(log2(numRow))/2 + avgWeightBit), param->validated); 
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculatePower((numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse));	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
				}
				// Read
				readDynamicEnergyArray = 0;
				readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs activityColWrite
				//20241031 update change the vdd to access voltage
				readDynamicEnergyArray += capRow2 * cell.accessVoltage * cell.accessVoltage; // Selected WL
				readDynamicEnergyArray *= numRow * activityRowRead * numColMuxed;

				readDynamicEnergy = 0;
				readDynamicEnergy += wlDecoder.readDynamicEnergy;
				readDynamicEnergy += wlNewDecoderDriver.readDynamicEnergy;
				readDynamicEnergy += wlDecoderDriver.readDynamicEnergy;
				readDynamicEnergy += ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += sarADC.readDynamicEnergy;
				
				readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;
				readDynamicEnergyAccum = adder.readDynamicEnergy + dff.readDynamicEnergy + shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergyOther = wlDecoder.readDynamicEnergy + wlNewDecoderDriver.readDynamicEnergy + wlDecoderDriver.readDynamicEnergy + ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);

				// Write					
				// writeDynamicEnergyArray = writeDynamicEnergyArray;
				// writeDynamicEnergy = 0;
				// if (cell.writeVoltage > 1.5) {
					// wllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// bllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// sllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);	
					// writeDynamicEnergy += wllevelshifter.writeDynamicEnergy + bllevelshifter.writeDynamicEnergy + sllevelshifter.writeDynamicEnergy; 
				// }
				// writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
				// writeDynamicEnergy += wlNewDecoderDriver.writeDynamicEnergy;
				// writeDynamicEnergy += wlDecoderDriver.writeDynamicEnergy;
				// writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;
				
				// Leakage
				leakage = 0;
				leakage += wlDecoder.leakage;
				leakage += wlDecoderDriver.leakage;
				leakage += wlNewDecoderDriver.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				leakage += dff.leakage;
				leakage += adder.leakage;
				leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
					
			cout << "\n================ RRAM/FeFET (Sequential) Block Energy & Power ================" << endl;
				 cout << "--- Read Dynamic Energy (J) ---" << endl;
				 cout << "WL Decoder Read Energy:       " << wlDecoder.readDynamicEnergy << " J" << endl;
				 if (cell.accessType == CMOS_access) {
				     cout << "WL Decoder Driver Read E:     " << wlNewDecoderDriver.readDynamicEnergy << " J" << endl;
				 } else {
				     cout << "WL Decoder Driver Read E:     " << wlDecoderDriver.readDynamicEnergy << " J" << endl;
				 }
				 if (numColMuxed > 1) {
				     cout << "MUX Decoder Read Energy:      " << muxDecoder.readDynamicEnergy << " J" << endl;
				     cout << "MUX Read Energy:              " << mux.readDynamicEnergy << " J" << endl;
				 }
				 if (param->SARADC) {
				     cout << "SAR ADC Read Energy:          " << sarADC.readDynamicEnergy << " J" << endl;
				 } else {
				     cout << "Multilevel SA Read Energy:    " << multilevelSenseAmp.readDynamicEnergy << " J" << endl;
				     if (avgWeightBit > 1) cout << "Multilevel SA Enc Read Energy:" << multilevelSAEncoder.readDynamicEnergy << " J" << endl;
				 }
				 cout << "Adder Read Energy:            " << adder.readDynamicEnergy << " J" << endl;
				 cout << "DFF Read Energy:              " << dff.readDynamicEnergy << " J" << endl;
				 if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Read Energy:  " << shiftAddWeight.readDynamicEnergy << " J" << endl;
				 if (numReadPulse > 1) cout << "ShiftAdd Input Read Energy:   " << shiftAddInput.readDynamicEnergy << " J" << endl;

				 cout << "--- Write Dynamic Energy (J) ---" << endl;
				 cout << "WL Decoder Write Energy:      " << wlDecoder.writeDynamicEnergy << " J" << endl;
				 if (cell.accessType == CMOS_access) {
				     cout << "WL Decoder Driver Write E:    " << wlNewDecoderDriver.writeDynamicEnergy << " J" << endl;
				 } else {
				     cout << "WL Decoder Driver Write E:    " << wlDecoderDriver.writeDynamicEnergy << " J" << endl;
				 }
				 cout << "SL Switch Matrix Write Energy:" << slSwitchMatrix.writeDynamicEnergy << " J" << endl;

				 if (cell.writeVoltage > 1.5) {
				     cout << "WL Level Shifter Write Energy:" << wllevelshifter.writeDynamicEnergy << " J" << endl;
				     cout << "BL Level Shifter Write Energy:" << bllevelshifter.writeDynamicEnergy << " J" << endl;
				     cout << "SL Level Shifter Write Energy:" << sllevelshifter.writeDynamicEnergy << " J" << endl;
				 }

				 cout << "--- Leakage Power (W) ---" << endl;
				 cout << "WL Decoder Leakage:           " << wlDecoder.leakage << " W" << endl;
				 if (cell.accessType == CMOS_access) {
				     cout << "WL Decoder Driver Leakage:    " << wlNewDecoderDriver.leakage << " W" << endl;
				 } else {
				     cout << "WL Decoder Driver Leakage:    " << wlDecoderDriver.leakage << " W" << endl;
				 }
				 cout << "SL Switch Matrix Leakage:     " << slSwitchMatrix.leakage << " W" << endl;
				 if (numColMuxed > 1) {
				     cout << "MUX Decoder Leakage:          " << muxDecoder.leakage << " W" << endl;
				     cout << "MUX Leakage:                  " << mux.leakage << " W" << endl;
				 }
				 if (param->SARADC) {
				     cout << "SAR ADC Leakage:              " << sarADC.leakage << " W" << endl;
				 } else {
				     cout << "Multilevel SA Leakage:        " << multilevelSenseAmp.leakage << " W" << endl;
				     if (avgWeightBit > 1) cout << "Multilevel SA Enc Leakage:    " << multilevelSAEncoder.leakage << " W" << endl;
				 }
				 cout << "Adder Leakage:                " << adder.leakage << " W" << endl;
				 cout << "DFF Leakage:                  " << dff.leakage << " W" << endl;
				 if (numCellPerSynapse > 1) cout << "ShiftAdd Weight Leakage:      " << shiftAddWeight.leakage << " W" << endl;
				 if (numReadPulse > 1) cout << "ShiftAdd Input Leakage:       " << shiftAddInput.leakage << " W" << endl;
				 cout << "======================================================================" << endl;
			} else if (conventionalParallel) {
				double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double capBL = lengthCol * 0.2e-15/1e-6;

				// 1.4 update: ADC update
				// param->reference_energy_peri = capRow2/param->numColSubArray * tech.vdd * tech.vdd * numRow;
				// for RRAM, the reference columns can always be turned on

				if (cell.accessType == CMOS_access) {
					wlNewSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
				} else {
					wlSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				}
				slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				// Anni update: numAdd
				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numAdd);
				} else {
					multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
					multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
				}
				if (numAdd > 1) {	
					adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)); 
					dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
				}
				if (numCellPerSynapse > 1) {
					shiftAddWeight.CalculatePower((numCellPerSynapse-1)*ceil(numColMuxed/numCellPerSynapse));	
				}
				if (numReadPulse > 1) {
					shiftAddInput.CalculatePower(ceil(numColMuxed/numCellPerSynapse));		
				}
				// Read
				readDynamicEnergyArray = 0;
				// Anni update:  * numAdd
				// readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells * numAdd; // Selected BLs activityColWrite -> no need already considered in multilevelsenseamp
				//20241031 update change the vdd to access voltage
				readDynamicEnergyArray += capRow2 * cell.accessVoltage * cell.accessVoltage * numRow * activityRowRead; // Selected WL
				// 1.4 update: buffer insertion
				readDynamicEnergyArray += (drivecapin + drivecapout) * cell.accessVoltage * cell.accessVoltage * param->buffernumber * 2 * numRow * activityRowRead;
				readDynamicEnergyArray *= numColMuxed;
				
				readDynamicEnergy = 0;
				readDynamicEnergy += wlNewSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				// Anni update: adder, dff, mux
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += sarADC.readDynamicEnergy;	
				
				readDynamicEnergyADC = readDynamicEnergyArray + multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy + sarADC.readDynamicEnergy;


				// Anni update: accum, other
				readDynamicEnergyAccum = shiftAddWeight.readDynamicEnergy + shiftAddInput.readDynamicEnergy + adder.readDynamicEnergy + dff.readDynamicEnergy;				
				readDynamicEnergyOther = wlNewSwitchMatrix.readDynamicEnergy + wlSwitchMatrix.readDynamicEnergy + ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				
				// Write				
				// writeDynamicEnergyArray = writeDynamicEnergyArray;
				// writeDynamicEnergy = 0;
				// if (cell.writeVoltage > 1.5) {
					// wllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// bllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// sllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// writeDynamicEnergy += wllevelshifter.writeDynamicEnergy + bllevelshifter.writeDynamicEnergy + sllevelshifter.writeDynamicEnergy;
				// }
				// writeDynamicEnergy += wlNewSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += wlSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage = 0;
				// 1.4 update: repeater leakage
				double repeater_leakage = CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * param->buffernumber *2 * param->numRowSubArray;
				leakage += repeater_leakage ;											   
				leakage += wlSwitchMatrix.leakage;
				leakage += wlNewSwitchMatrix.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				leakage += shiftAddWeight.leakage + shiftAddInput.leakage;
				// Anni update
				leakage += dff.leakage;
				leakage += adder.leakage;
				
			 cout << "\n================ RRAM/FeFET Block Energy & Power ================" << endl;
        		 cout << "--- Read Dynamic Energy (J) ---" << endl;
        		 if (wlDecoder.initialized) cout << "WL Decoder Read Energy:       " << wlDecoder.readDynamicEnergy << " J" << endl;
        		 if (wlNewDecoderDriver.initialized) cout << "WL New Decoder Driver Read E: " << wlNewDecoderDriver.readDynamicEnergy << " J" << endl;
        		 if (wlDecoderDriver.initialized) cout << "WL Decoder Driver Read Energy:" << wlDecoderDriver.readDynamicEnergy << " J" << endl;
        		 if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Read Energy: " << wlSwitchMatrix.readDynamicEnergy << " J" << endl;
        		 if (wlNewSwitchMatrix.initialized) cout << "WL New Switch Matrix Read E:  " << wlNewSwitchMatrix.readDynamicEnergy << " J" << endl;
        		 if (slSwitchMatrix.initialized) cout << "SL Switch Matrix Read Energy: " << slSwitchMatrix.readDynamicEnergy << " J" << endl;

        		 if (precharger.initialized) {
        		     cout << "Precharger Read Energy:       " << precharger.readDynamicEnergy << " J" << endl;
        		 }
        		 if (numColMuxed > 1) {
        		     if (muxDecoder.initialized) cout << "MUX Decoder Read Energy:      " << muxDecoder.readDynamicEnergy << " J" << endl;
        		     if (mux.initialized) cout << "MUX Read Energy:              " << mux.readDynamicEnergy << " J" << endl;
        		 }
        		 if (param->SARADC) {
        		     if (sarADC.initialized) cout << "SAR ADC Read Energy:          " << sarADC.readDynamicEnergy << " J" << endl;
        		 } else {
        		     if (multilevelSenseAmp.initialized) cout << "Multilevel SA Read Energy:    " << multilevelSenseAmp.readDynamicEnergy << " J" << endl;
        		 }
        		 if (numAdd > 1) {
        		     if (adder.initialized) cout << "Adder Read Energy:            " << adder.readDynamicEnergy << " J" << endl;
        		     if (dff.initialized) cout << "DFF Read Energy:              " << dff.readDynamicEnergy << " J" << endl;
        		 }
        		 if (numCellPerSynapse > 1) {
        		     if (shiftAddWeight.initialized) cout << "ShiftAdd Weight Read Energy:  " << shiftAddWeight.readDynamicEnergy << " J" << endl;
        		 }
        		 if (numReadPulse > 1) {
        		     if (shiftAddInput.initialized) cout << "ShiftAdd Input Read Energy:   " << shiftAddInput.readDynamicEnergy << " J" << endl;
        		 }

        		 cout << "--- Write Dynamic Energy (J) ---" << endl;
        		 if (wlDecoder.initialized) cout << "WL Decoder Write Energy:      " << wlDecoder.writeDynamicEnergy << " J" << endl;
        		 if (wlNewDecoderDriver.initialized) cout << "WL New Decoder Driver Write E:" << wlNewDecoderDriver.writeDynamicEnergy << " J" << endl;
        		 if (wlDecoderDriver.initialized) cout << "WL Decoder Driver Write E:    " << wlDecoderDriver.writeDynamicEnergy << " J" << endl;
        		 if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Write Energy:" << wlSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		 if (wlNewSwitchMatrix.initialized) cout << "WL New Switch Matrix Write E: " << wlNewSwitchMatrix.writeDynamicEnergy << " J" << endl;
        		 if (slSwitchMatrix.initialized) cout << "SL Switch Matrix Write Energy:" << slSwitchMatrix.writeDynamicEnergy << " J" << endl;

        		 if (cell.writeVoltage > 1.5) {
        		     if (wllevelshifter.initialized) cout << "WL Level Shifter Write Energy:" << wllevelshifter.writeDynamicEnergy << " J" << endl;
        		     if (bllevelshifter.initialized) cout << "BL Level Shifter Write Energy:" << bllevelshifter.writeDynamicEnergy << " J" << endl;
        		     if (sllevelshifter.initialized) cout << "SL Level Shifter Write Energy:" << sllevelshifter.writeDynamicEnergy << " J" << endl;
        		 }

        		 cout << "--- Leakage Power (W) ---" << endl;
        		 if (wlDecoder.initialized) cout << "WL Decoder Leakage:           " << wlDecoder.leakage << " W" << endl;
        		 if (wlNewDecoderDriver.initialized) cout << "WL New Decoder Drv Leakage:   " << wlNewDecoderDriver.leakage << " W" << endl;
        		 if (wlSwitchMatrix.initialized) cout << "WL Switch Matrix Leakage:     " << wlSwitchMatrix.leakage << " W" << endl;
        		 if (slSwitchMatrix.initialized) cout << "SL Switch Matrix Leakage:     " << slSwitchMatrix.leakage << " W" << endl;
        		 if (numColMuxed > 1 && muxDecoder.initialized) cout << "MUX Decoder Leakage:          " << muxDecoder.leakage << " W" << endl;
        		 if (param->SARADC && sarADC.initialized) cout << "SAR ADC Leakage:              " << sarADC.leakage << " W" << endl;
        		 cout << "=================================================================" << endl;
			
			} else if (BNNsequentialMode || XNORsequentialMode) {
				double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
				double numWriteCells = (int)ceil((double)numCol/*numWriteColMuxed*/); 
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double capBL = lengthCol * 0.2e-15/1e-6;
			
				wlDecoder.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				if (cell.accessType == CMOS_access) {
					wlNewDecoderDriver.CalculatePower(numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				} else {
					wlDecoderDriver.CalculatePower(numReadCells, numWriteCells, numRow*activityRowRead*numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite);
				}
				slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				
				// 1.4 update 230615
				// rowCurrentSenseAmp.CalculatePower(columnResistance, numRow*activityRowRead);
		
				multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
				multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);

				adder.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells);
				// Anni update
				dff.CalculatePower(numColMuxed*numRow*activityRowRead, numReadCells*(ceil(log2(numRow))/2+1), param->validated); 
				
				// Read
				readDynamicEnergyArray = 0;
				readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells; // Selected BLs activityColWrite
				//20241031 update change the vdd to access voltage
				readDynamicEnergyArray += capRow2 * cell.accessVoltage * cell.accessVoltage; // Selected WL
				readDynamicEnergyArray *= numRow * activityRowRead * numColMuxed;

				readDynamicEnergy = 0;
				readDynamicEnergy += wlDecoder.readDynamicEnergy;
				readDynamicEnergy += wlNewDecoderDriver.readDynamicEnergy;
				readDynamicEnergy += wlDecoderDriver.readDynamicEnergy;
				// Anni update
				readDynamicEnergy += ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy + multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;

				// Write				
				// writeDynamicEnergyArray = writeDynamicEnergyArray;
				// writeDynamicEnergy = 0;
				// if (cell.writeVoltage > 1.5) {
					// wllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// bllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// sllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// writeDynamicEnergy += wllevelshifter.writeDynamicEnergy + bllevelshifter.writeDynamicEnergy + sllevelshifter.writeDynamicEnergy;
				// }
				// writeDynamicEnergy += wlDecoder.writeDynamicEnergy;
				// writeDynamicEnergy += wlNewDecoderDriver.writeDynamicEnergy;
				// writeDynamicEnergy += wlDecoderDriver.writeDynamicEnergy;
				// writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage = 0;
				leakage += wlDecoder.leakage;
				leakage += wlDecoderDriver.leakage;
				leakage += wlNewDecoderDriver.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += rowCurrentSenseAmp.leakage;
				leakage += dff.leakage;
				leakage += adder.leakage;
				
			} else if (BNNparallelMode || XNORparallelMode) {
				double numReadCells = (int)ceil((double)numCol/numColMuxed);    // similar parameter as numReadCellPerOperationNeuro, which is for SRAM
				int numWriteOperationPerRow = (int)ceil((double)numCol*activityColWrite/numWriteCellPerOperationNeuro);
				double capBL = lengthCol * 0.2e-15/1e-6;

				if (cell.accessType == CMOS_access) {
					wlNewSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
				} else {
					wlSwitchMatrix.CalculatePower(numColMuxed, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				}
				slSwitchMatrix.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead, activityColWrite);
				if (numColMuxed > 1) {
					mux.CalculatePower(numColMuxed);	// Mux still consumes energy during row-by-row read
					muxDecoder.CalculatePower(numColMuxed, 1);
				}
				// Anni update: numAdd
				if (param->SARADC) {
					sarADC.CalculatePower(columnResistance, numAdd);
				} else {
					multilevelSenseAmp.CalculatePower(columnResistance, numAdd);
					multilevelSAEncoder.CalculatePower(numColMuxed * numAdd);
				}
				if (numAdd > 1) {	
					adder.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)); 
					dff.CalculatePower(numColMuxed*(numAdd-1), ceil(numCol/numColMuxed)*(log2(levelOutput) + ceil(log2(numAdd))/2), param->validated);
				}
				// Read
				readDynamicEnergyArray = 0;
				// Anni update

				
				// readDynamicEnergyArray += capBL * cell.readVoltage * cell.readVoltage * numReadCells * numAdd; // Selected BLs activityColWrite -> Already considered in multilevelsenseamp.cpp
				//20241031 update change the vdd to access voltage
				readDynamicEnergyArray += capRow2 * cell.accessVoltage * cell.accessVoltage * numRow * activityRowRead; // Selected WL
				readDynamicEnergyArray *= numColMuxed;
				// 1.4 update: ADC update
				// param->reference_energy_peri = capRow2/param->numColSubArray * tech.vdd * tech.vdd * numRow;

				readDynamicEnergy = 0;
				readDynamicEnergy += wlNewSwitchMatrix.readDynamicEnergy;
				readDynamicEnergy += wlSwitchMatrix.readDynamicEnergy;
				// Anni update
				readDynamicEnergy += ((numColMuxed > 1)==true? (mux.readDynamicEnergy + muxDecoder.readDynamicEnergy):0);
				readDynamicEnergy += adder.readDynamicEnergy;
				readDynamicEnergy += dff.readDynamicEnergy;
				readDynamicEnergy += multilevelSenseAmp.readDynamicEnergy;
				readDynamicEnergy += multilevelSAEncoder.readDynamicEnergy;
				readDynamicEnergy += readDynamicEnergyArray;
				readDynamicEnergy += sarADC.readDynamicEnergy;

				// Write				
				// writeDynamicEnergyArray = writeDynamicEnergyArray;
				// writeDynamicEnergy = 0;
				// if (cell.writeVoltage > 1.5) {
					// wllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// bllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// sllevelshifter.CalculatePower(0, 2*numWriteOperationPerRow*numRow*activityRowWrite, activityRowRead);
					// writeDynamicEnergy += wllevelshifter.writeDynamicEnergy + bllevelshifter.writeDynamicEnergy + sllevelshifter.writeDynamicEnergy;
				// }
				// writeDynamicEnergy += wlNewSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += wlSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += slSwitchMatrix.writeDynamicEnergy;
				// writeDynamicEnergy += writeDynamicEnergyArray;				
				
				// Leakage
				leakage = 0;
				leakage += wlSwitchMatrix.leakage;
				leakage += wlNewSwitchMatrix.leakage;
				leakage += slSwitchMatrix.leakage;
				leakage += ((numColMuxed > 1)==true? (mux.leakage):0);
				leakage += ((numColMuxed > 1)==true? (muxDecoder.leakage):0);
				leakage += multilevelSenseAmp.leakage;
				leakage += multilevelSAEncoder.leakage;
				// Anni update
				leakage += dff.leakage;
				leakage += adder.leakage;

			}
		} 
	}
}

void SubArray::PrintProperty() {

	if (cell.memCellType == Type::SRAM) {

		cout << endl << endl;
	    cout << "Array:" << endl;
	    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
	    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;

		precharger.PrintProperty("precharger");
		sramWriteDriver.PrintProperty("sramWriteDriver");

		if (conventionalSequential) {
			wlDecoder.PrintProperty("wlDecoder");
			senseAmp.PrintProperty("senseAmp");
			dff.PrintProperty("dff");
			adder.PrintProperty("adder");
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		} else if (conventionalParallel) {
			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		} else if (BNNsequentialMode || XNORsequentialMode) {
			wlDecoder.PrintProperty("wlDecoder");
			senseAmp.PrintProperty("senseAmp");
			dff.PrintProperty("dff");
			adder.PrintProperty("adder");
		} else if (BNNparallelMode || XNORparallelMode) {
			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
		} else {
			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		}

	} else if (cell.memCellType == Type::_2TnC) {
    		    cout << endl << endl;
    		    cout << "Array:" << endl;
    		    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
		    cout << "  [3D] staircase = " << areaStaircase*1e12 << " um^2,  core (array + staircase) = " << areaCore*1e12 << " um^2" << endl;
		    //cout << "  [3D] mode (0=CNA/1=CUA/2=CBA): " << param->integrationMode << endl;
		    //cout << "  [3D] staircase / via routing: " << areaStaircase*1e12 << " / " << areaViaRouting*1e12 << " um^2" << endl;
		    //cout << "  [3D] memory / logic die:      " << areaMemoryDie*1e12 << " / " << areaLogicDie*1e12 << " um^2" << endl;
		    //cout << "  [3D] bond pad / chip:         " << areaBondPad*1e12 << " / " << chipFootprint*1e12 << " um^2" << endl;
		    //cout << "  [3D] array efficiency:        " << (areaArray/chipFootprint)*100 << " %" << endl;
    		    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
    		    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;

    		    wwlSwitchMatrix.PrintProperty("wwlSwitchMatrix");
    		    wplSwitchMatrix.PrintProperty("wplSwitchMatrix");
    		    rblSwitchMatrix.PrintProperty("rblSwitchMatrix");
    		    wblSwitchMatrix.PrintProperty("wblSwitchMatrix");
    		    sslSwitchMatrix.PrintProperty("sslSwitchMatrix");
    		    rslSwitchMatrix.PrintProperty("rslSwitchMatrix");

    		    if (numColMuxed > 1) {
    		        mux.PrintProperty("mux");
    		        muxDecoder.PrintProperty("muxDecoder");
    		    }
    		    if (param->SARADC) {
    		        sarADC.PrintProperty("sarADC");
    		    } else {
    		        multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
    		        multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
    		    }
		    currentSenseAmp.PrintProperty("currentSenseAmp");
    		    if (numReadPulse > 1) {
    		        shiftAddWeight.PrintProperty("shiftAddWeight");
    		        shiftAddInput.PrintProperty("shiftAddInput");
    		    }

	} else if (cell.memCellType == Type::_1TnC) {
                  cout << endl << endl;
                  cout << "Array:" << endl;
                  cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
                  cout << "  [3D] staircase = " << areaStaircase*1e12 << " um^2,  core (array + staircase) = " << areaCore*1e12 << " um^2" << endl;
		  // cout << "  [3D] mode (0=CNA/1=CUA/2=CBA): " << param->integrationMode << endl;
                  // cout << "  [3D] staircase / via routing: " << areaStaircase*1e12 << " / " << areaViaRouting*1e12 << " um^2" << endl;
                  // cout << "  [3D] memory / logic die:      " << areaMemoryDie*1e12 << " / " << areaLogicDie*1e12 << " um^2" << endl;
                  // cout << "  [3D] bond pad / chip:         " << areaBondPad*1e12 << " / " << chipFootprint*1e12 << " um^2" << endl;
		  // cout << "  [3D] array efficiency:        " << (areaArray/chipFootprint)*100 << " %" << endl;
		  cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
                  cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;

                  wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
                  plSwitchMatrix.PrintProperty("plSwitchMatrix");
                  blSwitchMatrix.PrintProperty("blSwitchMatrix");

                  if (numColMuxed > 1) {
                      mux.PrintProperty("mux");
                      muxDecoder.PrintProperty("muxDecoder");
                  }
                  if (param->SARADC) {
                      sarADC.PrintProperty("sarADC");
                  } else {
                      multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
                      multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
                  }
		  currentSenseAmp.PrintProperty("currentSenseAmp");
                  if (numReadPulse > 1) {
                      shiftAddWeight.PrintProperty("shiftAddWeight");
                      shiftAddInput.PrintProperty("shiftAddInput");
                  }

	} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {

		cout << endl << endl;
	    	cout << "Array:" << endl;
	    	cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
	    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
		cout << "Write Latency = " << writeLatencyArray*1e9 << "ns" << endl;

		if (conventionalSequential) {
			wlDecoder.PrintProperty("wlDecoder");
			if (cell.accessType == CMOS_access) {
				wlNewDecoderDriver.PrintProperty("wlNewDecoderDriver");
			} else {
				wlDecoderDriver.PrintProperty("wlDecoderDriver");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp or single-bit SenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
			adder.PrintProperty("adder");
			dff.PrintProperty("dff");
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		} else if (conventionalParallel) {
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			} else {
				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		} else if (BNNsequentialMode || XNORsequentialMode) {
			wlDecoder.PrintProperty("wlDecoder");
			if (cell.accessType == CMOS_access) {
				wlNewDecoderDriver.PrintProperty("wlNewDecoderDriver");
			} else {
				wlDecoderDriver.PrintProperty("wlDecoderDriver");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			rowCurrentSenseAmp.PrintProperty("currentSenseAmp");
			adder.PrintProperty("adder");
			dff.PrintProperty("dff");
		} else if (BNNparallelMode || XNORparallelMode) {
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			} else {
				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
		} else {
			if (cell.accessType == CMOS_access) {
				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
			} else {
				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
			}
			slSwitchMatrix.PrintProperty("slSwitchMatrix");
			mux.PrintProperty("mux");
			muxDecoder.PrintProperty("muxDecoder");
			if (param->SARADC) {
				sarADC.PrintProperty("sarADC");
			} else {
				multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
				multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
			}
			if (numReadPulse > 1) {
				shiftAddWeight.PrintProperty("shiftAddWeight");
				shiftAddInput.PrintProperty("shiftAddInput");
			}
		}
	}
	FunctionUnit::PrintProperty("SubArray");
	cout << "Used Area = " << usedArea*1e12 << "um^2" << endl;
	cout << "Empty Area = " << emptyArea*1e12 << "um^2" << endl;
}

// void SubArray::PrintProperty() {
// 
// 	if (cell.memCellType == Type::SRAM) {
// 		
// 		cout << endl << endl;
// 	    cout << "Array:" << endl;
// 	    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
// 	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
// 	    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
// 		
// 		precharger.PrintProperty("precharger");
// 		sramWriteDriver.PrintProperty("sramWriteDriver");
// 		
// 		if (conventionalSequential) {
// 			wlDecoder.PrintProperty("wlDecoder");			
// 			senseAmp.PrintProperty("senseAmp");
// 			dff.PrintProperty("dff"); 
// 			adder.PrintProperty("adder");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		} else if (conventionalParallel) {
// 			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		} else if (BNNsequentialMode || XNORsequentialMode) {
// 			wlDecoder.PrintProperty("wlDecoder");			
// 			senseAmp.PrintProperty("senseAmp");
// 			dff.PrintProperty("dff"); 
// 			adder.PrintProperty("adder");
// 		} else if (BNNparallelMode || XNORparallelMode) {
// 			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 		} else {
// 			wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		}
// 		
// 	} else if (cell.memCellType == Type::_2TnC) {
//     		    cout << endl << endl;
//     		    cout << "Array:" << endl;
//     		    cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
//     		    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
//     		    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
//     		    
//     		    wwlSwitchMatrix.PrintProperty("wwlSwitchMatrix");
//     		    wplSwitchMatrix.PrintProperty("wplSwitchMatrix");
//     		    rblSwitchMatrix.PrintProperty("rblSwitchMatrix");
//     		    wblSwitchMatrix.PrintProperty("wblSwitchMatrix");
//     		    sslSwitchMatrix.PrintProperty("sslSwitchMatrix");
//     		    rslSwitchMatrix.PrintProperty("rslSwitchMatrix");
//     		    
//     		    if (numColMuxed > 1) {
//     		        mux.PrintProperty("mux");
//     		        muxDecoder.PrintProperty("muxDecoder");
//     		    }
//     		    if (param->SARADC) {
//     		        sarADC.PrintProperty("sarADC");
//     		    } else {
//     		        multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
//     		        multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
//     		    }
//     		    if (numReadPulse > 1) {
//     		        shiftAddWeight.PrintProperty("shiftAddWeight");
//     		        shiftAddInput.PrintProperty("shiftAddInput");
//     		    }
//     		
// 	} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
// 		
// 		cout << endl << endl;
// 	    	cout << "Array:" << endl;
// 	    	cout << "Area = " << heightArray*1e6 << "um x " << widthArray*1e6 << "um = " << areaArray*1e12 << "um^2" << endl;
// 	    cout << "Read Dynamic Energy = " << readDynamicEnergyArray*1e12 << "pJ" << endl;
// 	    cout << "Write Dynamic Energy = " << writeDynamicEnergyArray*1e12 << "pJ" << endl;
// 		cout << "Write Latency = " << writeLatencyArray*1e9 << "ns" << endl;
// 
// 		if (conventionalSequential) {
// 			wlDecoder.PrintProperty("wlDecoder");
// 			if (cell.accessType == CMOS_access) {
// 				wlNewDecoderDriver.PrintProperty("wlNewDecoderDriver");
// 			} else {
// 				wlDecoderDriver.PrintProperty("wlDecoderDriver");
// 			} 
// 			slSwitchMatrix.PrintProperty("slSwitchMatrix");
// 			mux.PrintProperty("mux");
// 			muxDecoder.PrintProperty("muxDecoder");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp or single-bit SenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 			adder.PrintProperty("adder");
// 			dff.PrintProperty("dff");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		} else if (conventionalParallel) {
// 			if (cell.accessType == CMOS_access) {
// 				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
// 			} else {
// 				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			}
// 			slSwitchMatrix.PrintProperty("slSwitchMatrix");
// 			mux.PrintProperty("mux");
// 			muxDecoder.PrintProperty("muxDecoder");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		} else if (BNNsequentialMode || XNORsequentialMode) {
// 			wlDecoder.PrintProperty("wlDecoder");
// 			if (cell.accessType == CMOS_access) {
// 				wlNewDecoderDriver.PrintProperty("wlNewDecoderDriver");
// 			} else {
// 				wlDecoderDriver.PrintProperty("wlDecoderDriver");
// 			} 
// 			slSwitchMatrix.PrintProperty("slSwitchMatrix");
// 			mux.PrintProperty("mux");
// 			muxDecoder.PrintProperty("muxDecoder");
// 			rowCurrentSenseAmp.PrintProperty("currentSenseAmp");
// 			adder.PrintProperty("adder");
// 			dff.PrintProperty("dff");
// 		} else if (BNNparallelMode || XNORparallelMode) {
// 			if (cell.accessType == CMOS_access) {
// 				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
// 			} else {
// 				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			}
// 			slSwitchMatrix.PrintProperty("slSwitchMatrix");
// 			mux.PrintProperty("mux");
// 			muxDecoder.PrintProperty("muxDecoder");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 		} else {
// 			if (cell.accessType == CMOS_access) {
// 				wlNewSwitchMatrix.PrintProperty("wlNewSwitchMatrix");
// 			} else {
// 				wlSwitchMatrix.PrintProperty("wlSwitchMatrix");
// 			}
// 			slSwitchMatrix.PrintProperty("slSwitchMatrix");
// 			mux.PrintProperty("mux");
// 			muxDecoder.PrintProperty("muxDecoder");
// 			multilevelSenseAmp.PrintProperty("multilevelSenseAmp");
// 			multilevelSAEncoder.PrintProperty("multilevelSAEncoder");
// 			if (numReadPulse > 1) {
// 				shiftAddWeight.PrintProperty("shiftAddWeight");
// 				shiftAddInput.PrintProperty("shiftAddInput");
// 			}
// 		}
// 	} 
// 	FunctionUnit::PrintProperty("SubArray");
// 	cout << "Used Area = " << usedArea*1e12 << "um^2" << endl;
// 	cout << "Empty Area = " << emptyArea*1e12 << "um^2" << endl;

