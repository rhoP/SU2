/*!
 * \file CDiscAdjTurbMLSolver.hpp
 * \brief Headers of the CDiscAdjSolver class
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

#pragma once

#include "CSolver.hpp"
#include "../variables/CDiscAdjVariable.hpp"

/*!
 * \class CDiscAdjTurbMLSolver
 * \brief Main class for defining the discrete adjoint solver for turbulence modeling with ML.
 * \ingroup Discrete_Adjoint
 * \author R. Pochampalli
 */
class CDiscAdjTurbMLSolver final : public CSolver {
private:
    unsigned short KindDirect_Solver;
    CSolver *direct_solver;
    su2double ObjFunc_Value;       /*!< \brief Value of the objective function. */
    su2double Reg_Value;           /*!< \brief Value of the regularization function. */
    su2double Reg_Param_Value;     /*!< \brief Value of the objective function. */

    vector<su2double> Turb_Params;
    su2double *Sensitivity_Turb_params = nullptr; /*!< \brief Auxiliary vector for the geometry solution (dimension nDim instead of nVar). */

    CDiscAdjVariable* nodes = nullptr;  /*!< \brief The highest level in the variable hierarchy this solver can safely use. */

    /*!
     * \brief Return nodes to allow CSolver::base_nodes to be set.
     */
    inline CVariable* GetBaseClassPointerToNodes() override { return nodes; }

public:

    /*!
     * \overload
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] solver - Initialize the discrete adjoint solver with the corresponding direct solver.
     * \param[in] Kind_Solver - The kind of direct solver.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
     */
    CDiscAdjTurbMLSolver(CGeometry *geometry, CConfig *config, CSolver* solver, unsigned short Kind_Solver, unsigned short iMesh);

    /*!
     * \brief Destructor of the class.
     */
    ~CDiscAdjTurbMLSolver(void);

    /*!
     * \brief Performs the preprocessing of the adjoint AD-based solver.
     *        Registers all necessary variables on the tape. Called while tape is active.
     * \param[in] geometry_container - The geometry container holding all grid levels.
     * \param[in] config_container - The particular config.
     */
    void RegisterSolution(CGeometry *geometry, CConfig *config) override;

    /*!
     * \brief Performs the preprocessing of the adjoint AD-based solver.
     *        Registers all necessary variables that are output variables on the tape.
     *        Called while tape is active.
     * \param[in] geometry_container - The geometry container holding all grid levels.
     * \param[in] config_container - The particular config.
     */
    void RegisterOutput(CGeometry *geometry, CConfig *config) override;

    /*!
     * \brief Sets the adjoint values of the output of the flow (+turb.) iteration
     *         before evaluation of the tape.
     * \param[in] geometry - The geometrical definition of the problem.
     * \param[in] config - The particular config.
     */
    void SetAdjoint_Output(CGeometry *geometry, CConfig *config) override;

    /*!
     * \brief Sets the adjoint values of the input variables of the flow (+turb.) iteration
     *        after tape has been evaluated.
     * \param[in] geometry - The geometrical definition of the problem.
     * \param[in] config - The particular config.
     */
    void ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config) override;



    /*!
     * \brief Register the objective function as output.
     * \param[in] geometry - The geometrical definition of the problem.
     */
    void RegisterObj_Func(CConfig *config) override;

    /*!
     * \brief Set the objective function.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void SetAdj_ObjFunc(CGeometry *geometry, CConfig* config) override;
    /*!
     * \brief Prepare the solver for a new recording.
     * \param[in] kind_recording - Kind of AD recording.
     */
    void SetRecording(CGeometry *geometry, CConfig *config) override;

    /*!
     * \brief A virtual member.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     * \param[in] reset - If true reset variables to their initial values.
     */
    void RegisterVariables(CGeometry *geometry,
                           CConfig *config,
                           bool reset = false) override;

    /*!
     * \brief A virtual member.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    void ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config) override;

    /*!
     * \brief Update the dual-time derivatives.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container vector with all the solutions.
     * \param[in] config - Definition of the particular problem.
     * \param[in] iMesh - Index of the mesh in multigrid computations.
     * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
     * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
     * \param[in] Output - boolean to determine whether to print output.
     */
    void Preprocessing(CGeometry *geometry,
                       CSolver **solver_container,
                       CConfig *config,
                       unsigned short iMesh,
                       unsigned short iRKStep,
                       unsigned short RunTime_EqSystem,
                       bool Output) override;

    /*!
     * \brief Load a solution from a restart file.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver - Container vector with all of the solvers.
     * \param[in] config - Definition of the particular problem.
     * \param[in] val_iter - Current external iteration number.
     * \param[in] val_update_geo - Flag for updating coords and grid velocity.
     */
    void LoadRestart(CGeometry **geometry,
                     CSolver ***solver,
                     CConfig *config,
                     int val_iter,
                     bool val_update_geo) override;

    /*!
 * \brief Value of the objective function
 * \param[in] config - Definition of the particular problem.
 */
   virtual su2double Get_Objective_Value(CConfig *config);

    /*!
       * \brief Get parameter sensitivity.
       * \param[in] point_index: index of the point.
       * \param[out] returns the sensitivity of the indexed ML parameter.
       */
    virtual su2double GetMLParamSens(unsigned long point_index) override {return Sensitivity_Turb_params[point_index];}
};
