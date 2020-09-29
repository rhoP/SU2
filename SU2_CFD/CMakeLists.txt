

project(SU2_CFD
LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 14)




include(${CMAKE_SOURCE_DIR}/cmake/functions.cmake)

set(CMAKE_MODULE_PATH "/home/rohit/libtorch/share/cmake/Torch;${CMAKE_MODULE_PATH}")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")



set(SU2_CORE_SOURCES
./src/CMarkerProfileReaderFVM.cpp
./src/SU2_CFD.cpp
./src/definition_structure.cpp
./src/drivers/CDiscAdjMultizoneDriver.cpp
./src/drivers/CDriver.cpp
./src/drivers/CDummyDriver.cpp
./src/drivers/CMultizoneDriver.cpp
./src/drivers/CSinglezoneDriver.cpp
./src/drivers/CDiscAdjSinglezoneDriver.cpp
./src/fluid_model.cpp
./src/fluid_model_inc.cpp
./src/fluid_model_pig.cpp
./src/fluid_model_ppr.cpp
./src/fluid_model_pvdw.cpp
./src/integration/CFEM_DG_Integration.cpp
./src/integration/CIntegration.cpp
./src/integration/CIntegrationFactory.cpp
./src/integration/CMultiGridIntegration.cpp
./src/integration/CSingleGridIntegration.cpp
./src/integration/CStructuralIntegration.cpp
./src/interfaces/CInterface.cpp
./src/interfaces/cfd/CConservativeVarsInterface.cpp
./src/interfaces/cfd/CMixingPlaneInterface.cpp
./src/interfaces/cfd/CSlidingInterface.cpp
./src/interfaces/cht/CConjugateHeatInterface.cpp
./src/interfaces/fsi/CDiscAdjFlowTractionInterface.cpp
./src/interfaces/fsi/CDisplacementsInterface.cpp
./src/interfaces/fsi/CFlowTractionInterface.cpp
./src/limiters/CLimiterDetails.cpp
./src/numerics/CNumerics.cpp
./src/numerics/continuous_adjoint/adj_convection.cpp
./src/numerics/continuous_adjoint/adj_diffusion.cpp
./src/numerics/continuous_adjoint/adj_sources.cpp
./src/numerics/elasticity/CFEAElasticity.cpp
./src/numerics/elasticity/CFEALinearElasticity.cpp
./src/numerics/elasticity/CFEANonlinearElasticity.cpp
./src/numerics/elasticity/nonlinear_models.cpp
./src/numerics/flow/convection/ausm_slau.cpp
./src/numerics/flow/convection/centered.cpp
./src/numerics/flow/convection/cusp.cpp
./src/numerics/flow/convection/fds.cpp
./src/numerics/flow/convection/fvs.cpp
./src/numerics/flow/convection/hllc.cpp
./src/numerics/flow/convection/roe.cpp
./src/numerics/flow/flow_diffusion.cpp
./src/numerics/flow/flow_sources.cpp
./src/numerics/heat.cpp
./src/numerics/radiation.cpp
./src/numerics/template.cpp
./src/numerics/transition.cpp
./src/numerics/turbulent/turb_convection.cpp
./src/numerics/turbulent/turb_diffusion.cpp
./src/numerics/turbulent/turb_sources.cpp
./src/output/CAdjElasticityOutput.cpp
./src/output/CAdjFlowCompOutput.cpp
./src/output/CAdjFlowIncOutput.cpp
./src/output/CAdjHeatOutput.cpp
./src/output/CBaselineOutput.cpp
./src/output/CElasticityOutput.cpp
./src/output/CFlowCompFEMOutput.cpp
./src/output/CFlowIncOutput.cpp
./src/output/CFlowOutput.cpp
./src/output/CHeatOutput.cpp
./src/output/CMeshOutput.cpp
./src/output/CMultizoneOutput.cpp
./src/output/COutputFactory.cpp
./src/output/filewriter/CCSVFileWriter.cpp
./src/output/filewriter/CFEMDataSorter.cpp
./src/output/filewriter/CFVMDataSorter.cpp
./src/output/filewriter/CParallelDataSorter.cpp
./src/output/filewriter/CParallelFileWriter.cpp
./src/output/filewriter/CParaviewBinaryFileWriter.cpp
./src/output/filewriter/CParaviewFileWriter.cpp
./src/output/filewriter/CParaviewVTMFileWriter.cpp
./src/output/filewriter/CParaviewXMLFileWriter.cpp
./src/output/filewriter/CSTLFileWriter.cpp
./src/output/filewriter/CSU2BinaryFileWriter.cpp
./src/output/filewriter/CSU2FileWriter.cpp
./src/output/filewriter/CSU2MeshFileWriter.cpp
./src/output/filewriter/CSurfaceFEMDataSorter.cpp
./src/output/filewriter/CSurfaceFVMDataSorter.cpp
./src/output/filewriter/CTecplotBinaryFileWriter.cpp
./src/output/filewriter/CTecplotFileWriter.cpp
./src/output/output_physics.cpp
./src/output/output_structure_legacy.cpp
./src/output/tools/CWindowingTools.cpp
./src/output/COutput.cpp
./src/output/CFlowCompOutput.cpp
./src/python_wrapper_structure.cpp
./src/solvers/CAdjEulerSolver.cpp
./src/solvers/CAdjNSSolver.cpp
./src/solvers/CAdjTurbSolver.cpp
./src/solvers/CBaselineSolver.cpp
./src/solvers/CBaselineSolver_FEM.cpp
./src/solvers/CDiscAdjFEASolver.cpp
./src/solvers/CDiscAdjMeshSolver.cpp
./src/solvers/CFEASolver.cpp
./src/solvers/CFEM_DG_EulerSolver.cpp
./src/solvers/CFEM_DG_NSSolver.cpp
./src/solvers/CHeatSolver.cpp
./src/solvers/CIncEulerSolver.cpp
./src/solvers/CIncNSSolver.cpp
./src/solvers/CMeshSolver.cpp
./src/solvers/CNSSolver.cpp
./src/solvers/CRadP1Solver.cpp
./src/solvers/CRadSolver.cpp
./src/solvers/CSolver.cpp
./src/solvers/CSolverFactory.cpp
./src/solvers/CTemplateSolver.cpp
./src/solvers/CTransLMSolver.cpp
./src/solvers/CTurbSSTSolver.cpp
./src/solvers/CTurbSolver.cpp
./src/solvers/CEulerSolver.cpp
./src/solvers/CDiscAdjSolver.cpp
./src/solvers/CTurbSASolver.cpp
./src/transport_model.cpp
./src/variables/CAdjEulerVariable.cpp
./src/variables/CAdjNSVariable.cpp
./src/variables/CAdjTurbVariable.cpp
./src/variables/CBaselineVariable.cpp
./src/variables/CDiscAdjFEABoundVariable.cpp
./src/variables/CDiscAdjFEAVariable.cpp
./src/variables/CDiscAdjMeshBoundVariable.cpp
./src/variables/CEulerVariable.cpp
./src/variables/CFEABoundVariable.cpp
./src/variables/CFEAVariable.cpp
./src/variables/CHeatVariable.cpp
./src/variables/CIncEulerVariable.cpp
./src/variables/CIncNSVariable.cpp
./src/variables/CMeshBoundVariable.cpp
./src/variables/CMeshElement.cpp
./src/variables/CMeshVariable.cpp
./src/variables/CNSVariable.cpp
./src/variables/CRadP1Variable.cpp
./src/variables/CRadVariable.cpp
./src/variables/CTransLMVariable.cpp
./src/variables/CTurbSSTVariable.cpp
./src/variables/CVariable.cpp
./src/variables/CTurbSAVariable.cpp
./src/variables/CDiscAdjVariable.cpp
./src/variables/CTurbVariable.cpp
./src/iteration_structure.cpp)


set(SU2_CFD_SOURCES include/SU2_CFD.hpp src/SU2_CFD.cpp )


set(SU2_CORE_LIB SU2Core)
set(SU2_CFD_EXE SU2_CFD)



add_library(${SU2_CORE_LIB} STATIC ${SU2_CORE_SOURCES})
target_include_directories(${SU2_CORE_LIB} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/include)
target_include_directories(${SU2_CORE_LIB} PUBLIC ~/libtorch/include)

add_library(SU2::SU2Core ALIAS ${SU2_CORE_LIB})

add_executable(${SU2_CFD_EXE} ${SU2_CFD_SOURCES})

set_property(TARGET ${SU2_CFD_EXE} PROPERTY CXX_STANDARD 14)


target_link_libraries(${SU2_CFD_EXE} PRIVATE "${TORCH_LIBRARIES}")


su2_setup_target(${SU2_CORE_LIB})
su2_setup_target(${SU2_CFD_EXE} LIBS SU2::SU2Core SU2::SU2 WITH_THREADS)
