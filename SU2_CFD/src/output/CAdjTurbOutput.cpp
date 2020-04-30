/*!
 * \file CAdjTurbOutput.cpp
 * \brief Main subroutines of the adjoint turbulence modeling with machine learning output class.
 * \author R. Pochampalli
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */


#include "../../include/output/CAdjTurbOutput.hpp"

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"

CAdjTurbOutput::CAdjTurbOutput(CConfig *config, unsigned short nDim) : COutput(config, nDim, false) {


    /*--- Set the default history fields if nothing is set in the config file ---*/

    if (nRequestedHistoryFields == 0){
        requestedHistoryFields.emplace_back("ITER");
        requestedHistoryFields.emplace_back("RMS_RES");
        requestedHistoryFields.emplace_back("SENSITIVITY");
        requestedHistoryFields.emplace_back("OBJ_FUNC");
        nRequestedHistoryFields = requestedHistoryFields.size();
    }

    if (nRequestedScreenFields == 0){
        if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
        requestedScreenFields.emplace_back("INNER_ITER");
        requestedScreenFields.emplace_back("RMS_ADJ_PRESSURE");
        requestedScreenFields.emplace_back("RMS_ADJ_NU_TILDE");
        requestedScreenFields.emplace_back("SENS_FIELD");
        requestedHistoryFields.emplace_back("OBJ_FUNC");
        nRequestedScreenFields = requestedScreenFields.size();
    }

    if (nRequestedVolumeFields == 0){
        requestedVolumeFields.emplace_back("COORDINATES");
        requestedVolumeFields.emplace_back("SOLUTION");
        requestedVolumeFields.emplace_back("SENSITIVITY");
        nRequestedVolumeFields = requestedVolumeFields.size();
    }

    stringstream ss;
    ss << "Zone " << config->GetiZone() << " (Adj. turbulence)";
    multiZoneHeaderString = ss.str();

    /*--- Set the volume filename --- */

    volumeFilename = config->GetAdj_FileName();

    /*--- Set the surface filename --- */

    surfaceFilename = config->GetSurfAdjCoeff_FileName();

    /*--- Set the restart filename --- */

    restartFilename = config->GetRestart_AdjFileName();

    /*--- Add the obj. function extension --- */

    restartFilename = config->GetObjFunc_Extension(restartFilename);

    /*--- Add the field sensitivity output filename ---*/

    fieldSensitivityFileName = config->Get_FieldSensitivity_FileName();

    /*--- Set the default convergence field --- */

    if (convFields.empty() ) convFields.emplace_back("RMS_ADJ_PRESSURE");

    /*--- Set the field sensitivity file name if modeling turbulence ---*/
    FieldSensitivityFileName = config->Get_FieldSensitivity_FileName();

}


void CAdjTurbOutput::SetHistoryOutputFields(CConfig *config){

    //Log Residuals
    AddHistoryOutput("RMS_ADJ_PRESSURE",    "rms[A_P]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint Pressure.", HistoryFieldType::RESIDUAL);

    AddHistoryOutput("RMS_ADJ_NU_TILDE", "rms[A_nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);

    //Maximum Residuals
    AddHistoryOutput("MAX_ADJ_PRESSURE",    "max[A_Rho]",  ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint Pressure.", HistoryFieldType::RESIDUAL);

    AddHistoryOutput("MAX_ADJ_NU_TILDE", "max[A_nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);

    //Sensitivities
    AddHistoryOutput("SENS_FIELD", "Sens[Beta]", ScreenOutputFormat::SCIENTIFIC, "SENSITIVITY", "Total Sensitivity of the field parameters.", HistoryFieldType::COEFFICIENT);

    //Objective Function
    AddHistoryOutput("OBJ_FUNC", "ObjFun", ScreenOutputFormat::SCIENTIFIC, "OBJ_FUNC", "", HistoryFieldType::COEFFICIENT);

}

inline void CAdjTurbOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver) {

    SetHistoryOutputValue("RMS_ADJ_PRESSURE", log10(solver[ADJFLOW_SOL]->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_ADJ_NU_TILDE", log10(solver[ADJTURB_SOL]->GetRes_RMS(0)));

    SetHistoryOutputValue("MAX_ADJ_NU_TILDE", log10(solver[ADJTURB_SOL]->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_ADJ_PRESSURE", log10(solver[ADJFLOW_SOL]->GetRes_Max(0)));

    SetHistoryOutputValue("SENS_FIELD", (solver[ADJTURB_SOL]->GetTotalFieldSens()));

    SetHistoryOutputValue("OBJ_FUNC", solver[FLOW_SOL]->GetTotal_ComboObj());

}

void CAdjTurbOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){

    CVariable* Node_AdjTurb = solver[ADJTURB_SOL]->GetNodes();
    CVariable* Node_AdjFlow = solver[ADJFLOW_SOL]->GetNodes();

    CPoint*    Node_Geo  = geometry->node[iPoint];


    SetVolumeOutputValue("COORD-X", iPoint,  Node_Geo->GetCoord(0));
    SetVolumeOutputValue("COORD-Y", iPoint,  Node_Geo->GetCoord(1));
    if (nDim == 3)
        SetVolumeOutputValue("COORD-Z", iPoint, Node_Geo->GetCoord(2));

    SetVolumeOutputValue("ADJ_PRESSURE",   iPoint, Node_AdjFlow->GetSolution(iPoint, 0));
    SetVolumeOutputValue("RES_ADJ_PRESSURE",   iPoint, Node_AdjFlow->GetSolution(iPoint, 0) - Node_AdjFlow->GetSolution_Old(iPoint, 0));

    SetVolumeOutputValue("ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
    SetVolumeOutputValue("RES_ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
    SetVolumeOutputValue("BETA", iPoint, *(solver[ADJTURB_SOL]->Get_iParamML(iPoint)));

    SetVolumeOutputValue("ADJ_BETA", iPoint, solver[ADJTURB_SOL]->GetMLParamSens(iPoint));

}

void CAdjTurbOutput::SetVolumeOutputFields(CConfig *config){
    // Grid coordinates
    AddVolumeOutput("COORD-X", "x", "COORDINATES", "x-component of the coordinate vector");
    AddVolumeOutput("COORD-Y", "y", "COORDINATES", "y-component of the coordinate vector");
    if (nDim == 3)
        AddVolumeOutput("COORD-Z", "z", "COORDINATES", "z-component of the coordinate vector");

    /// BEGIN_GROUP: SOLUTION, DESCRIPTION: The SOLUTION variables of the adjoint solver.
    /// DESCRIPTION: Adjoint Pressure.
    AddVolumeOutput("ADJ_PRESSURE",    "Adjoint_Pressure",    "SOLUTION", "Adjoint pressure");

    /// BEGIN_GROUP: RESIDUAL, DESCRIPTION: Residuals of the SOLUTION variables.
    /// DESCRIPTION: Residual of the adjoint Pressure.
    AddVolumeOutput("RES_ADJ_PRESSURE",    "Residual_Adjoint_Pressure",    "RESIDUAL", "Residual of the adjoint pressure");


    AddVolumeOutput("ADJ_NU_TILDE", "Adjoint_Nu_Tilde", "SOLUTION", "Adjoint Spalart-Allmaras variable");
    AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL", "Residual of the adjoint Spalart-Allmaras variable");
    AddVolumeOutput("BETA", "Field_Parameters", "SOLUTION", "field parameter");

    AddVolumeOutput("ADJ_BETA", "Field_Sensitivity", "SENSITIVITY", "sensitivity of the field parameter");

}

bool CAdjTurbOutput::SetInit_Residuals(CConfig *config){

    return (config->GetTime_Marching() != STEADY && (curInnerIter == 0))||
           (config->GetTime_Marching() == STEADY && (curTimeIter < 2));

}


bool CAdjTurbOutput::SetUpdate_Averages(CConfig *config){
    return false;

//  return (config->GetUnsteady_Simulation() != STEADY && !dualtime);

}

void CAdjTurbOutput::LoadSurfaceData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex){


}