/*!
 * \file CTurbSASolver.hpp
 * \brief Headers of the CTurbSASolver class
 * \author A. Bueno.
 * \version 7.0.4 "Blackbird"
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

#include "CTurbSolver.hpp"
#include <torch/script.h>


class PicElem{
private:
    su2double x;
    su2double y;
    su2double total_kernel;

public:
    vector<unsigned long> neighbors;
    vector<unsigned long> kernels;
    PicElem(){
        x = 0.;
        y = 0.;
    }
    ~PicElem()=default;
    PicElem(su2double x_val, su2double y_val){
        x = x_val;
        y = y_val;
    }
    su2double get_x() const{
        return x;
    }
    su2double get_y() const{
        return y;
    }
    void set_x(su2double x_val){
        x = x_val;
    }
    void set_y(su2double y_val){
        y = y_val;
    }
    void translate(su2double xt, su2double yt){
        x -= xt;
        y -= yt;
    }
    void set_total_kernel(su2double ValKernel){
        total_kernel = ValKernel;
    }
    su2double get_total_kernel(){
        return total_kernel;
    }
};



/*!
 * \class CTurbSASolver
 * \brief Main class for defining the turbulence model solver.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSASolver final : public CTurbSolver {
private:
  su2double nu_tilde_Inf, nu_tilde_Engine, nu_tilde_ActDisk;

  su2double field_param_DD;

  vector<vector<unsigned long>> neighbors;

  vector<vector<vector<PicElem>>> picture_kernels;

  vector<unsigned long> domain_t;
  su2double kernel_parameter{1.0E-3};
  su2double nbDistance{0.03};
  su2double nbRadius{0.02};


  torch::jit::script::Module module;

  vector<vector<PicElem>> baseCoords;


  /*!
   * \brief A virtual member.
   * \param[in] solver - Solver container
   * \param[in] geometry - Geometrical definition.
   * \param[in] config - Definition of the particular problem.
   */
  void SetDES_LengthScale(CSolver** solver,
                          CGeometry *geometry,
                          CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   */
  CTurbSASolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] FluidModel
   */
  CTurbSASolver(CGeometry *geometry, CConfig *config, unsigned short iMesh, CFluidModel* FluidModel);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSASolver(void) override;

  /*!
   * \brief Restart residual and compute gradients.
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
   * \brief A virtual member.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   */
  void Postprocessing(CGeometry *geometry,
                      CSolver **solver_container,
                      CConfig *config,
                      unsigned short iMesh) override;
  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics_container - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics **numerics_container,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Template(CGeometry *geometry,
                       CSolver **solver_container,
                       CNumerics *numerics,
                       CConfig *config,
                       unsigned short iMesh) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry *geometry,
                          CSolver **solver_container,
                          CNumerics *conv_numerics,
                          CNumerics *visc_numerics,
                          CConfig *config,
                          unsigned short val_marker) override;

  /*!
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry *geometry,
                    CSolver **solver_container,
                    CNumerics *conv_numerics,
                    CNumerics *visc_numerics,
                    CConfig *config,
                    unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry *geometry,
                CSolver **solver_container,
                CNumerics *conv_numerics,
                CNumerics *visc_numerics,
                CConfig *config,
                unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_Turbo(CGeometry *geometry,
                      CSolver **solver_container,
                      CNumerics *conv_numerics,
                      CNumerics *visc_numerics,
                      CConfig *config,
                      unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet_MixingPlane(CGeometry *geometry,
                            CSolver **solver_container,
                            CNumerics *conv_numerics,
                            CNumerics *visc_numerics,
                            CConfig *config,
                            unsigned short val_marker) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry *geometry,
                 CSolver **solver_container,
                 CNumerics *conv_numerics,
                 CNumerics *visc_numerics,
                 CConfig *config,
                 unsigned short val_marker) override;

  /*!
   * \brief Impose the engine inflow boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Inflow(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose the engine exhaust boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Engine_Exhaust(CGeometry *geometry,
                         CSolver **solver_container,
                         CNumerics *conv_numerics,
                         CNumerics *visc_numerics,
                         CConfig *config,
                         unsigned short val_marker) override;

  /*!
   * \brief Impose the interface boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Interface_Boundary(CGeometry *geometry,
                             CSolver **solver_container,
                             CNumerics *numerics,
                             CConfig *config,
                             unsigned short val_marker) override;

  /*!
   * \brief Impose the fluid interface boundary condition using tranfer data.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void BC_Fluid_Interface(CGeometry *geometry,
                          CSolver **solver_container,
                          CNumerics *conv_numerics,
                          CNumerics *visc_numerics,
                          CConfig *config) override;

  /*!
   * \brief Impose the near-field boundary condition using the residual.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_NearField_Boundary(CGeometry *geometry,
                             CSolver **solver_container,
                             CNumerics *numerics,
                             CConfig *config,
                             unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Inlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_ActDisk_Outlet(CGeometry *geometry,
                        CSolver **solver_container,
                        CNumerics *conv_numerics,
                        CNumerics *visc_numerics,
                        CConfig *config,
                        unsigned short val_marker) override;

  /*!
   * \brief Impose an actuator disk inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] val_inlet_surface - Boolean for whether val_marker is an inlet
   */
  void BC_ActDisk(CGeometry *geometry,
                  CSolver **solver_container,
                  CNumerics *conv_numerics,
                  CNumerics *visc_numerics,
                  CConfig *config,
                  unsigned short val_marker,
                  bool val_inlet_surface) override;

  /*!
   * \brief Set the solution using the Freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  inline void SetFreeStream_Solution(CConfig *config) override {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) nodes->SetSolution(iPoint, 0, nu_tilde_Inf);
  }

  /*!
   * \brief Store of a set of provided inlet profile values at a vertex.
   * \param[in] val_inlet - vector containing the inlet values for the current vertex.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   * \param[in] iVertex - Vertex of the marker <i>iMarker</i> where the inlet is being set.
   */
  void SetInletAtVertex(su2double *val_inlet,
                        unsigned short iMarker,
                        unsigned long iVertex) override;

  /*!
   * \brief Get the set of value imposed at an inlet.
   * \param[in] val_inlet - vector returning the inlet values for the current vertex.
   * \param[in] val_inlet_point - Node index where the inlet is being set.
   * \param[in] val_kind_marker - Enumerated type for the particular inlet type.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param config - Definition of the particular problem.
   * \return Value of the face area at the vertex.
   */
  su2double GetInletAtVertex(su2double *val_inlet,
                             unsigned long val_inlet_point,
                             unsigned short val_kind_marker,
                             string val_marker,
                             CGeometry *geometry,
                             CConfig *config) const override;

  /*!
   * \brief Set a uniform inlet profile
   *
   * The values at the inlet are set to match the values specified for
   * inlets in the configuration file.
   *
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMarker - Surface marker where the coefficient is computed.
   */
  void SetUniformInlet(CConfig* config, unsigned short iMarker) override;

  /*!
   * \brief Get the value of nu tilde at the far-field.
   * \return Value of nu tilde at the far-field.
   */
  inline su2double GetNuTilde_Inf(void) const override { return nu_tilde_Inf; }

  /*!
   * \brief Compute nu tilde from the wall functions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void SetNuTilde_WF(CGeometry *geometry,
                     CSolver **solver_container,
                     CNumerics *conv_numerics,
                     CNumerics *visc_numerics,
                     CConfig *config,
                     unsigned short val_marker) override;

  void ReadFieldParameters(CConfig *config, CGeometry *geometry);


  su2double GetFieldRegularization(CConfig *config, CGeometry *geometry) override;

  /*
   * \brief Computes the kernels for each point and stores them into the kernel vector
   * \param[in] geometry - Geometrical definition of the problem
   * \param[in] config - Definition of the particular problem
   */
  void SetKernels(CConfig *config, CGeometry *geometry);

    /*
     * \brief returns the kernel vector for a point
     * \param[in] iPoint - index of the point
     */
  // vector<su2double> GetKernelValue(unsigned long iPoint){ return kernels[iPoint];}

  /*
   * \brief returns the kernel parameter
   * \param[out] the distance parameter in the kernel
   */
  // su2double GetKernelParameter(){return kernel_parameter;}


  /*
   * \brief Sets the kernel parameter to some value
   * \param[in] val_parameter - distance to set the kernel parameter to
   */
/*  void SetKernelParameter(su2double val_parameter){
      kernel_parameter = val_parameter;
  }*/

/*
  void SetNbDistanceKernelReg(su2double dist){
      nbDistance = dist;
  }
*/

/*  su2double GetNbDistanceKernelReg(){
      return nbDistance;
  }*/



  void SetNeighbors(CConfig *config, CGeometry *geometry);


  void SetTurbulenceModelCorrectionDomain(CConfig *config, CGeometry *geometry);

  void GenerateChannels(torch::Tensor& channels,  unsigned long iPoint, CSolver** solver,
          CNumerics* numerics, CGeometry* geometry);

  // vector<PicElem> GetNeighbors(unsigned long iPoint);

  su2double euclidean_distance(CGeometry* geometry, unsigned long iPoint, unsigned long jPoint){
        return pow(pow(geometry->node[iPoint]->GetCoord(0) - geometry->node[jPoint]->GetCoord(0), 2) +
                   pow(geometry->node[iPoint]->GetCoord(1) - geometry->node[jPoint]->GetCoord(1), 2), 0.5);
  }

  su2double n12(su2double xc) {
      auto t1 = 0.298222773 * sqrt(xc) - 0.127125232 * xc - 0.357907906 * pow(xc, 2) +
                0.291984971 * pow(xc, 3) - 0.105174606 * pow(xc, 4);
      return 0.594689181 * t1;
  }

  su2double n21(su2double xc) {
      return 1.05 * (0.2969 * sqrt(xc) - 0.1260 * xc - 0.3516 * pow(xc, 2) + 0.2843 * pow(xc, 3) - 0.1015 * pow(xc, 4));
  }


  void Compute_BaseCoordinates();

  vector<vector<su2double>> channel_stats = {{0.0, 0.6055165138445859}, {0.001868454135732041, 13207.48470087414},
                                             {2.753034850436209e-09, 42714.921478879245}, {2.753034850436209e-09, 37897.135874986765},
                                             {-12758.01543289031, 9518.125427962812}, {-16815.03497503484, 13932.639725076951},
                                             {-42714.921478879245, 1006.599410577314}, {-9520.534595526595, 12742.92357026478},
                                             {-3814.9494511811267, 657.7961930513623}, {-4805.409451433434, 1050.8959221561079},
                                             {-4565.9310762351615, 238.6522048116658}, {-983.8827420516634, 1903.761364052603},
                                             {-2253.838004301427, 9603.70130055727}, {0.0, 0.0012092927100644}};



};
