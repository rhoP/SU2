#include "../../include/solvers/CSolver.hpp"
#include "../../include/solvers/CSolverFactory.hpp"
#include "../../include/solvers/CEulerSolver.hpp"
#include "../../include/solvers/CIncEulerSolver.hpp"
#include "../../include/solvers/CNSSolver.hpp"
#include "../../include/solvers/CIncNSSolver.hpp"
#include "../../include/solvers/CTurbSASolver.hpp"
#include "../../include/solvers/CTurbSSTSolver.hpp"
#include "../../include/solvers/CTransLMSolver.hpp"
#include "../../include/solvers/CAdjEulerSolver.hpp"
#include "../../include/solvers/CAdjNSSolver.hpp"
#include "../../include/solvers/CAdjTurbSolver.hpp"
#include "../../include/solvers/CHeatSolverFVM.hpp"
#include "../../include/solvers/CFEASolver.hpp"
#include "../../include/solvers/CTemplateSolver.hpp"
#include "../../include/solvers/CDiscAdjSolver.hpp"
#include "../../include/solvers/CDiscAdjFEASolver.hpp"
#include "../../include/solvers/CFEM_DG_EulerSolver.hpp"
#include "../../include/solvers/CFEM_DG_NSSolver.hpp"
#include "../../include/solvers/CMeshSolver.hpp"
#include "../../include/solvers/CDiscAdjMeshSolver.hpp"

CSolver** CSolverFactory::createSolverContainer(ENUM_SOLVER kindSolver, CConfig *config, CGeometry *geometry, int iMGLevel){
  
  CSolver** solver;
  
  solver = new CSolver*[MAX_SOLS];
  for (unsigned int iSol = 0; iSol < MAX_SOLS; iSol++)
    solver[iSol] = nullptr;
  
  const bool allocDirect  = false;
  const bool allocAdjoint = true;
  
  ENUM_TURB_MODEL kindTurbModel = static_cast<ENUM_TURB_MODEL>(config->GetKind_Turb_Model());

  switch (kindSolver) {
    case TEMPLATE_SOLVER:
      solver[FLOW_SOL] = new CTemplateSolver(geometry, config);
      break;
    case INC_EULER:
      solver[FLOW_SOL] = new CIncEulerSolver(geometry, config, iMGLevel);
      break;
    case EULER:
      solver[FLOW_SOL] = new CEulerSolver(geometry, config, iMGLevel);
      break;
    case INC_NAVIER_STOKES:
      solver[FLOW_SOL] = new CIncNSSolver(geometry, config, iMGLevel);
      solver[HEAT_SOL] = createHeatSolver(solver, geometry,config, iMGLevel, allocDirect);
      break;
    case NAVIER_STOKES:
      solver[FLOW_SOL] = new CNSSolver(geometry, config, iMGLevel);
      break;
    case RANS:
      solver[FLOW_SOL] = new CNSSolver(geometry, config, iMGLevel);
      solver[TURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      break;
    case INC_RANS:
      solver[FLOW_SOL] = new CIncNSSolver(geometry, config, iMGLevel);
      solver[HEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect);
      solver[TURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      break;
    case HEAT_EQUATION_FVM:
      solver[HEAT_SOL] = new CHeatSolverFVM(geometry, config, iMGLevel);
      break;
    case ADJ_EULER:
      solver[FLOW_SOL]    = new CEulerSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CAdjEulerSolver(geometry, config, iMGLevel);
      break;
    case ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = new CNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CAdjNSSolver(geometry, config, iMGLevel);
      break;    
    case ADJ_RANS:
      solver[FLOW_SOL]    = new CNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CAdjNSSolver(geometry, config, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_EULER:
      solver[FLOW_SOL]    = new CEulerSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case DISC_ADJ_NAVIER_STOKES:
      solver[FLOW_SOL]    = new CNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case DISC_ADJ_RANS:
      solver[FLOW_SOL]    = new CNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_INC_EULER:
      solver[FLOW_SOL]    = new CIncEulerSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case DISC_ADJ_INC_NAVIER_STOKES:
      solver[FLOW_SOL]       = new CIncNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL]    = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      solver[HEAT_SOL]       = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJHEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_INC_RANS:
      solver[FLOW_SOL]    = new CIncNSSolver(geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      solver[TURB_SOL]    = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJTURB_SOL] = createTurbSolver(kindTurbModel, solver, geometry, config, iMGLevel, allocAdjoint);
      solver[HEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJHEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case DISC_ADJ_HEAT:
      solver[HEAT_SOL]    = createHeatSolver(solver, geometry, config, iMGLevel, allocDirect);
      solver[ADJHEAT_SOL] = createHeatSolver(solver, geometry, config, iMGLevel, allocAdjoint);
      break;
    case FEM_ELASTICITY:
      solver[FEA_SOL] = new CFEASolver(geometry, config);
      break;
    case DISC_ADJ_FEM:
      solver[FEA_SOL]    = new CFEASolver(geometry, config);
      solver[ADJFEA_SOL] = new CDiscAdjFEASolver(geometry, config, solver[FEA_SOL], RUNTIME_FEA_SYS, iMGLevel);
      break;
    case FEM_EULER:
      solver[FLOW_SOL] = createDGSolver(FEM_EULER, geometry, config, iMGLevel);
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      solver[FLOW_SOL] = createDGSolver(FEM_NAVIER_STOKES, geometry, config, iMGLevel);
      break;
    case FEM_RANS:
      solver[FLOW_SOL] = createDGSolver(FEM_RANS, geometry, config, iMGLevel);
      break;
    case DISC_ADJ_FEM_EULER:
      solver[FLOW_SOL]    = createDGSolver(FEM_EULER, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case DISC_ADJ_FEM_NS:
      solver[FLOW_SOL]    = createDGSolver(FEM_NAVIER_STOKES, geometry, config, iMGLevel);
      solver[ADJFLOW_SOL] = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_FLOW_SYS, iMGLevel);
      break;
    case DISC_ADJ_FEM_RANS:
      SU2_MPI::Error("Finite element turbulence model not yet implemented.", CURRENT_FUNCTION);
      break;
     default:
      solver = nullptr;
  }
  
  solver[MESH_SOL]    = createMeshSolver(solver, geometry, config, iMGLevel, allocDirect);
  solver[ADJMESH_SOL] = createMeshSolver(solver, geometry, config, iMGLevel, allocAdjoint);
  
  return solver;
  
}

CSolver* CSolverFactory::createTurbSolver(ENUM_TURB_MODEL kindTurbModel, CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, int adjoint){
  
  CSolver *turbSolver = nullptr;
  
  if (!adjoint){
    switch (kindTurbModel) {
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        turbSolver = new CTurbSASolver(geometry, config, iMGLevel, solver[FLOW_SOL]->GetFluidModel());
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        break;
      case SST: case SST_SUST:
        turbSolver = new CTurbSSTSolver(geometry, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        turbSolver->Postprocessing(geometry, solver, config, iMGLevel);
        solver[FLOW_SOL]->Preprocessing(geometry, solver, config, iMGLevel, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
        break;
      default:
        SU2_MPI::Error("Unknown turbulence model", CURRENT_FUNCTION);
        break;
    }
  } else {
    
    if (config->GetDiscrete_Adjoint()){
      if (!config->GetFrozen_Visc_Disc())
        turbSolver = new CDiscAdjSolver(geometry, config, solver[FLOW_SOL], RUNTIME_TURB_SYS, iMGLevel);
    } else {
      if (!config->GetFrozen_Visc_Cont())
        turbSolver = new CAdjTurbSolver(geometry, config, iMGLevel);
    }
  }

  return turbSolver;
}

CSolver* CSolverFactory::createHeatSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){
  
  CSolver *heatSolver = nullptr;
  
  if (config->GetWeakly_Coupled_Heat()){
    if (adjoint){
      if (config->GetDiscrete_Adjoint()){
        heatSolver = new CDiscAdjSolver(geometry, config, solver[HEAT_SOL], RUNTIME_HEAT_SYS, iMGLevel);
      } else {
        SU2_MPI::Error("No continuous adjoint heat solver available.", CURRENT_FUNCTION);
      }
    } else {
      heatSolver = new CHeatSolverFVM(geometry, config, iMGLevel);
    }
  }
  return heatSolver;
  
}

CSolver* CSolverFactory::createMeshSolver(CSolver **solver, CGeometry *geometry, CConfig *config, int iMGLevel, bool adjoint){
  
  CSolver *meshSolver = nullptr;
  
  if (config->GetDeform_Mesh() && iMGLevel == MESH_0){
    meshSolver = new CMeshSolver(geometry, config);
    if (adjoint && config->GetDiscrete_Adjoint()){
      meshSolver = new CDiscAdjMeshSolver(geometry, config, solver[MESH_SOL]);
    }
  }
  return meshSolver;
  
}

CSolver* CSolverFactory::createDGSolver(ENUM_SOLVER kindDGSolver, CGeometry *geometry, CConfig *config, int iMGLevel){
  
  CSolver *DGSolver = nullptr;
  
  switch (kindDGSolver) {
    case FEM_EULER:
      if (config->GetKind_FEM_DG_Shock() == PERSSON){
        DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      } else {
        DGSolver = new CFEM_DG_EulerSolver(geometry, config, iMGLevel);
      }
      break;
    case FEM_NAVIER_STOKES: case FEM_LES:
      DGSolver = new CFEM_DG_NSSolver(geometry, config, iMGLevel);
      break;
    case FEM_RANS:
      SU2_MPI::Error("DG RANS solver not yet implemented.", CURRENT_FUNCTION);
      break;
    default:
      SU2_MPI::Error("Requested DG solver not found", CURRENT_FUNCTION);
      break;
  }
  
  return DGSolver;
}
