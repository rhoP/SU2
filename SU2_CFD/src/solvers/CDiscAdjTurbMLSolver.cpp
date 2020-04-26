/*!
 * \file CDiscAdjTurbMLSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint problem .
 * \author T. Albring, R. Pochampalli
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

#include "../../include/solvers/CDiscAdjTurbMLSolver.hpp"


CDiscAdjTurbMLSolver::CDiscAdjTurbMLSolver(CGeometry *geometry,
                                           CConfig *config, CSolver *direct_solver,
                                           unsigned short Kind_Solver, unsigned short iMesh)  : CSolver() {

    unsigned short iVar, iMarker, iDim;
    unsigned long iVertex;
    string text_line, mesh_filename;
    ifstream restart_file;
    string filename, AdjExt;

    adjoint = true;

    nVar = direct_solver->GetnVar();
    nDim = geometry->GetnDim();

    /*-- Store some information about direct solver ---*/
    this->KindDirect_Solver = Kind_Solver;
    this->direct_solver = direct_solver;


    nMarker      = config->GetnMarker_All();
    nPoint       = geometry->GetnPoint();
    nPointDomain = geometry->GetnPointDomain();

    /*--- Define some auxiliary vectors related to the residual ---*/

    Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
    Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
    Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

    /*--- Instantiate the turbulence parameter adjoint variables---*/
    Sensitivity_Turb_params.reserve(nPoint);

    Turb_Params.reserve(nPointDomain);


    for (unsigned long iPoint=0; iPoint < nPoint; iPoint++){
        Sensitivity_Turb_params.emplace_back(0.0);
        Turb_Params.emplace_back(0.0);
    }


    /*--- Define some auxiliary vectors related to the residual for problems with a BGS strategy---*/

    if (config->GetMultizone_Residual()){

        Residual_BGS      = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_BGS[iVar]      = 1.0;
        Residual_Max_BGS  = new su2double[nVar];     for (iVar = 0; iVar < nVar; iVar++) Residual_Max_BGS[iVar]  = 1.0;

        /*--- Define some structures for locating max residuals ---*/

        Point_Max_BGS       = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max_BGS[iVar] = 0;
        Point_Max_Coord_BGS = new su2double*[nVar];
        for (iVar = 0; iVar < nVar; iVar++) {
            Point_Max_Coord_BGS[iVar] = new su2double[nDim];
            for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord_BGS[iVar][iDim] = 0.0;
        }

    }


    /*--- Define some structures for locating max residuals ---*/

    Point_Max = new unsigned long[nVar] ();

    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
        Point_Max_Coord[iVar] = new su2double[nDim]();
    }

    /*--- Define some auxiliary vectors related to the solution ---*/

    Solution = new su2double[nVar];

    for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 1e-16;



    /*--- Initialize the discrete adjoint solution to zero everywhere. ---*/

    nodes = new CDiscAdjVariable(Solution, nPoint, nDim, nVar, config);
    SetBaseClassPointerToNodes();

    /*--- Set which points are vertices and allocate boundary data. ---*/

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
            iVertex = geometry->node[iPoint]->GetVertex(iMarker);
            if (iVertex >= 0) {
                nodes->Set_isVertex(iPoint,true);
                break;
            }
        }

    /*--- Store the direct solution ---*/

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
        nodes->SetSolution_Direct(iPoint, direct_solver->GetNodes()->GetSolution(iPoint));
    }


    switch(KindDirect_Solver){
        case RUNTIME_FLOW_SYS:
            SolverName = "ADJ.FLOW";
            break;
        case RUNTIME_HEAT_SYS:
            SolverName = "ADJ.HEAT";
            break;
        case RUNTIME_TURB_SYS:
            SolverName = "ADJ.TURB";
            break;
        case RUNTIME_RADIATION_SYS:
            SolverName = "ADJ.RAD";
            break;
        default:
            SolverName = "ADJ.SOL";
            break;
    }
}

CDiscAdjTurbMLSolver::~CDiscAdjTurbMLSolver(void) {

    unsigned short iMarker;



    if (nodes != nullptr) delete nodes;

}

void CDiscAdjTurbMLSolver::SetRecording(CGeometry* geometry, CConfig *config){

    bool time_n1_needed = config->GetTime_Marching() == DT_STEPPING_2ND;
    bool time_n_needed = (config->GetTime_Marching() == DT_STEPPING_1ST) || time_n1_needed;

    unsigned long iPoint;
    unsigned short iVar;

    /*--- Reset the solution to the initial (converged) solution ---*/

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        direct_solver->GetNodes()->SetSolution(iPoint, nodes->GetSolution_Direct(iPoint));
    }

    if (time_n_needed) {
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
            for (iVar = 0; iVar < nVar; iVar++) {
                AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
            }
        }
    }
    if (time_n1_needed) {
        for (iPoint = 0; iPoint < nPoint; iPoint++) {
            for (iVar = 0; iVar < nVar; iVar++) {
                AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n1(iPoint)[iVar]);
            }
        }
    }

    /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
     * when running the discrete adjoint solver. ---*/

    direct_solver->Jacobian.SetValZero();

    /*--- Set indices to zero ---*/

    RegisterVariables(geometry, config, true);

}


void CDiscAdjTurbMLSolver::RegisterSolution(CGeometry *geometry, CConfig *config) {

    bool time_n1_needed = (config->GetTime_Marching() == DT_STEPPING_2ND);
    bool time_n_needed  = (config->GetTime_Marching() == DT_STEPPING_1ST) || time_n1_needed;
    bool input          = true;
    bool push_index     = !config->GetMultizone_Problem();

    /*--- Register solution at all necessary time instances and other variables on the tape ---*/

    direct_solver->GetNodes()->RegisterSolution(input, push_index);

    if (time_n_needed)
        direct_solver->GetNodes()->RegisterSolution_time_n();

    if (time_n1_needed)
        direct_solver->GetNodes()->RegisterSolution_time_n1();
}

void CDiscAdjTurbMLSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset) {

    unsigned long iPoint, global_index;
    nPoint   =  geometry->GetnPoint();

    if (KindDirect_Solver == RUNTIME_TURB_SYS){
        for (iPoint = 0; iPoint < nPoint; iPoint++){
            global_index = geometry->node[iPoint]->GetGlobalIndex();
            Turb_Params[global_index] = *geometry->MLParam_Container->Get_iParamML(iPoint);
            AD::RegisterInput(Turb_Params[global_index]);
        }
    }

    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CDiscAdjTurbMLSolver::RegisterOutput(CGeometry *geometry, CConfig *config) {

    bool input        = false;
    bool push_index   = !config->GetMultizone_Problem();

    /*--- Register variables as output of the solver iteration ---*/

    direct_solver->GetNodes()->RegisterSolution(input, push_index);
}

void CDiscAdjTurbMLSolver::RegisterObj_Func(CConfig *config) {

    /*--- Here we can add new (scalar) objective functions ---*/
    switch (config->GetKind_ObjFunc()) {
            case DRAG_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CD();
                if (config->GetFixed_CL_Mode()) ObjFunc_Value -= config->GetdCD_dCL() * direct_solver->GetTotal_CL();
                if (config->GetFixed_CM_Mode()) ObjFunc_Value -= config->GetdCD_dCMy() * direct_solver->GetTotal_CMy();
                break;
            case LIFT_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CL();
                break;
            case SIDEFORCE_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CSF();
                break;
            case EFFICIENCY:
                ObjFunc_Value = direct_solver->GetTotal_CEff();
                break;
            case MOMENT_X_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CMx();
                break;
            case MOMENT_Y_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CMy();
                break;
            case MOMENT_Z_COEFFICIENT:
                ObjFunc_Value = direct_solver->GetTotal_CMz();
                break;
            case EQUIVALENT_AREA:
                ObjFunc_Value = direct_solver->GetTotal_CEquivArea();
                break;
            case BUFFET_SENSOR:
                ObjFunc_Value = direct_solver->GetTotal_Buffet_Metric();
                break;
            case TOTAL_HEATFLUX:
                ObjFunc_Value = direct_solver->GetTotal_HeatFlux();
                break;
            case INVERSE_DESIGN_ML:
                ObjFunc_Value = Get_Objective_Value(config);
                break;
            default:
                ObjFunc_Value = 0.0;
                break;
    }

        /*--- Template for new objective functions where TemplateObjFunction()
         *  is the routine that returns the obj. function value. The computation
         * must be done while the tape is active, i.e. between AD::StartRecording() and
         * AD::StopRecording() in DiscAdjMeanFlowIteration::Iterate(). The best place is somewhere
         * inside MeanFlowIteration::Iterate().
         *
         * case TEMPLATE_OBJECTIVE:
         *    ObjFunc_Value = TemplateObjFunction();
         *    break;
         * ---*/

    if (rank == MASTER_NODE) {
        AD::RegisterOutput(ObjFunc_Value);
        cout<< "The Obj Func Val is " << ObjFunc_Value << endl;
    }
}

void CDiscAdjTurbMLSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config) {

    bool time_stepping = config->GetTime_Marching() != STEADY;
    unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
    unsigned long TimeIter = config->GetTimeIter();
    su2double seeding = 1.0;

    if (time_stepping) {
        if (TimeIter < IterAvg_Obj) {
            seeding = 1.0/((su2double)IterAvg_Obj);
        }
        else {
            seeding = 0.0;
        }
    }

    if (rank == MASTER_NODE) {
        SU2_TYPE::SetDerivative(ObjFunc_Value, SU2_TYPE::GetValue(seeding));
    } else {
        SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
    }
}

void CDiscAdjTurbMLSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

    bool time_n1_needed = config->GetTime_Marching() == DT_STEPPING_2ND;
    bool time_n_needed = (config->GetTime_Marching() == DT_STEPPING_1ST) || time_n1_needed;
    bool multizone = config->GetMultizone_Problem();

    unsigned short iVar;
    unsigned long iPoint;
    su2double residual;

    /*--- Set Residuals to zero ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar,0.0);
        SetRes_Max(iVar,0.0,0);
    }

    /*--- Set the old solution ---*/

    if(!multizone) nodes->Set_OldSolution();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {

        /*--- Extract the adjoint solution ---*/

        if(config->GetMultizone_Problem()) {
            direct_solver->GetNodes()->GetAdjointSolution_LocalIndex(iPoint,Solution);
        }
        else {
            direct_solver->GetNodes()->GetAdjointSolution(iPoint,Solution);
        }

        /*--- Store the adjoint solution ---*/

        nodes->SetSolution(iPoint,Solution);
    }

    if (time_n_needed) {
        for (iPoint = 0; iPoint < nPoint; iPoint++) {

            /*--- Extract the adjoint solution at time n ---*/

            direct_solver->GetNodes()->GetAdjointSolution_time_n(iPoint,Solution);

            /*--- Store the adjoint solution at time n ---*/

            nodes->Set_Solution_time_n(iPoint,Solution);
        }
    }
    if (time_n1_needed) {
        for (iPoint = 0; iPoint < nPoint; iPoint++) {

            /*--- Extract the adjoint solution at time n-1 ---*/

            direct_solver->GetNodes()->GetAdjointSolution_time_n1(iPoint,Solution);

            /*--- Store the adjoint solution at time n-1 ---*/

            nodes->Set_Solution_time_n1(iPoint,Solution);
        }
    }

    /*--- Set the residuals ---*/

    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
            residual = nodes->GetSolution(iPoint,iVar) - nodes->GetSolution_Old(iPoint,iVar);

            AddRes_RMS(iVar,residual*residual);
            AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
        }
    }

    if(multizone) nodes->Set_OldSolution();

    SetResidual_RMS(geometry, config);

}

void CDiscAdjTurbMLSolver::ExtractAdjoint_Variables(CGeometry *geometry,
                                                    CConfig *config) {
    nPoint = geometry->GetnPoint();
    unsigned long global_index;

    if (KindDirect_Solver == RUNTIME_TURB_SYS){
        for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
            global_index = geometry->node[iPoint]->GetGlobalIndex();
            Sensitivity_Turb_params[global_index] = SU2_TYPE::GetDerivative(Turb_Params[global_index]);
        }
    }
}

void CDiscAdjTurbMLSolver::SetAdjoint_Output(CGeometry *geometry,
                                             CConfig *config) {

    bool dual_time = (config->GetTime_Marching() == DT_STEPPING_1ST ||
                      config->GetTime_Marching() == DT_STEPPING_2ND);
    bool fsi = config->GetFSI_Simulation();

    unsigned short iVar;
    unsigned long iPoint;

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] = nodes->GetSolution(iPoint,iVar);
        }
        if (dual_time) {
            for (iVar = 0; iVar < nVar; iVar++) {
                Solution[iVar] += nodes->GetDual_Time_Derivative(iPoint,iVar);
            }
        }
        direct_solver->GetNodes()->SetAdjointSolution(iPoint,Solution);
    }
}

void CDiscAdjTurbMLSolver::Preprocessing(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config_container,
                                         unsigned short iMesh,
                                         unsigned short iRKStep,
                                         unsigned short RunTime_EqSystem,
                                         bool Output) {
    bool dual_time_1st = (config_container->GetTime_Marching() == DT_STEPPING_1ST);
    bool dual_time_2nd = (config_container->GetTime_Marching() == DT_STEPPING_2ND);
    bool dual_time = (dual_time_1st || dual_time_2nd);
    su2double *solution_n, *solution_n1;
    unsigned long iPoint;
    unsigned short iVar;
    if (dual_time) {
        for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++) {
            solution_n = nodes->GetSolution_time_n(iPoint);
            solution_n1 = nodes->GetSolution_time_n1(iPoint);
            for (iVar=0; iVar < nVar; iVar++) {
                nodes->SetDual_Time_Derivative(iPoint, iVar, solution_n[iVar]+nodes->GetDual_Time_Derivative_n(iPoint, iVar));
                nodes->SetDual_Time_Derivative_n(iPoint,iVar, solution_n1[iVar]);
            }
        }
    }
}

void CDiscAdjTurbMLSolver::LoadRestart(CGeometry **geometry,
                                       CSolver ***solver,
                                       CConfig *config,
                                       int val_iter,
                                       bool val_update_geo) {

    unsigned short iVar, iMesh;
    unsigned long iPoint, index, iChildren, Point_Fine, counter;
    su2double Area_Children, Area_Parent, *Solution_Fine;
    string restart_filename, filename;

    bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
    bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
    bool rans = ((config->GetKind_Solver() == DISC_ADJ_RANS) || (config->GetKind_Solver() == DISC_ADJ_INC_RANS)) ;

    /*--- Restart the solution from file information ---*/

    filename = config->GetSolution_AdjFileName();
    restart_filename = config->GetObjFunc_Extension(filename);

    restart_filename = config->GetFilename(restart_filename, "", val_iter);


    /*--- Read and store the restart metadata. ---*/

//  Read_SU2_Restart_Metadata(geometry[MESH_0], config, true, restart_filename);

    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
        Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
        Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Read all lines in the restart file ---*/

    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- Skip coordinates ---*/
    unsigned short skipVars = geometry[MESH_0]->GetnDim();

    /*--- Skip flow adjoint variables ---*/
    if (KindDirect_Solver== RUNTIME_TURB_SYS) {
        if (compressible) {
            skipVars += nDim + 2;
        }
        if (incompressible) {
            skipVars += nDim + 2;
        }
    }

    /*--- Skip flow adjoint and turbulent variables ---*/
    if (KindDirect_Solver == RUNTIME_RADIATION_SYS) {
        if (compressible) skipVars += nDim + 2;
        if (incompressible) skipVars += nDim + 2;
        if (rans) skipVars += solver[MESH_0][TURB_SOL]->GetnVar();
    }

    /*--- Load data from the restart into correct containers. ---*/

    counter = 0;
    for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

        /*--- Retrieve local index. If this node from the restart file lives
         on the current processor, we will load and instantiate the vars. ---*/

        iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

        if (iPoint_Local > -1) {

            /*--- We need to store this point's data, so jump to the correct
             offset in the buffer of data from the restart file and load it. ---*/

            index = counter*Restart_Vars[1] + skipVars;
            for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
            nodes->SetSolution(iPoint_Local,Solution);
            iPoint_Global_Local++;

            /*--- Increment the overall counter for how many points have been loaded. ---*/
            counter++;
        }

    }

    /*--- Detect a wrong solution file ---*/

    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

    if (rbuf_NotMatching != 0) {
        SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                       string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Communicate the loaded solution on the fine grid before we transfer
     it down to the coarse levels. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
        for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
            Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
            for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
            for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
                Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
                Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
                Solution_Fine = solver[iMesh-1][ADJFLOW_SOL]->GetNodes()->GetSolution(Point_Fine);
                for (iVar = 0; iVar < nVar; iVar++) {
                    Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
                }
            }
            solver[iMesh][ADJFLOW_SOL]->GetNodes()->SetSolution(iPoint, Solution);
        }
    }

    /*--- Delete the class memory that is used to load the restart. ---*/

    if (Restart_Vars != nullptr) delete [] Restart_Vars;
    if (Restart_Data != nullptr) delete [] Restart_Data;
    Restart_Vars = nullptr; Restart_Data = nullptr;

}

su2double CDiscAdjTurbMLSolver::Get_Objective_Value(CConfig *config) {

return 0.5 * pow((direct_solver->GetTotal_CL() - 1.074902), 2);
}

/*!
 * \brief Get maximum parameter sensitivity.
 * \param[out] returns the maximum sensitivity of the ML parameters.
 */
su2double CDiscAdjTurbMLSolver::GetTotalFieldSens () {
    return accumulate(Sensitivity_Turb_params.begin(), Sensitivity_Turb_params.end(), decltype(Sensitivity_Turb_params)::value_type(0.0));
}

void CDiscAdjTurbMLSolver::SetSensitivity(CGeometry *geometry, CSolver **solver, CConfig *config) {
    ExtractAdjoint_Variables(geometry,config);
}
