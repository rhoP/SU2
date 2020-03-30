/*!
 * \file CTurbSAVariable.hpp
 * \brief Declaration of the variables of turbulence models with machine learning.
 * \author R. Pochampalli
 * \version 7.0.2 "Blackbird"
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

#include "CTurbVariable.hpp"


/*!
 * \class CTurbSA_MLVariable
 * \brief Class for defining the variables of the turbulence model with machine learning
 * \ingroup Turbulence_Model
 * \author Rohit P.
 */

class CTurbSA_MLVariable : public CTurbVariable {

private:
    su2double gamma_BC;             /*!< \brief Value of the intermittency for the BC trans. model. */
    su2double DES_LengthScale;
    su2double Vortex_Tilting;
    su2double ML_Param;             /*!< \brief Value of the machine learning parameter for turbulence modeling. */
    unsigned long original_index;   /*!< \brief Original index of node before RCM reordering */
public:
    /*!
     * \brief Constructor of the class.
     */
    CTurbSA_MLVariable(void);

    /*!
     * \overload
     * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
     * \param[in] val_muT  - The eddy viscosity
     * \param[in] val_nDim - Number of dimensions of the problem.
     * \param[in] val_nvar - Number of variables of the problem.
     * \param[in] config - Definition of the particular problem.
     */
    CTurbSA_MLVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim,
                       unsigned short val_nvar, CConfig *config, su2double val_mlparam, unsigned long val_index);

    /*!
     * \brief Destructor of the class.
     */
    ~CTurbSA_MLVariable(void);

    /*!
     * \brief Set the harmonic balance source term.
     * \param[in] val_var - Index of the variable.
     * \param[in] val_source - Value of the harmonic balance source term. for the index <i>val_var</i>.
     */
    inline void SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) {HB_Source[val_var] = val_source; }

    /*!
     * \brief Get the harmonic balance source term.
     * \param[in] val_var - Index of the variable.
     * \return Value of the harmonic balance source term for the index <i>val_var</i>.
     */
    inline su2double GetHarmonicBalance_Source(unsigned short val_var) {return HB_Source[val_var]; }

    /*!
     * \brief Get the intermittency of the BC transition model.
     * \return Value of the intermittency of the BC transition model.
     */
    inline su2double GetGammaBC(void) {return gamma_BC; }

    /*!
     * \brief Set the intermittency of the BC transition model.
     * \param[in] val_gamma - New value of the intermittency.
     */
    inline void SetGammaBC(su2double val_gamma) {gamma_BC = val_gamma; }

    /*!
     * \brief Get the DES length scale
     * \return Value of the DES length Scale.
     */
    inline su2double GetDES_LengthScale(void) {return DES_LengthScale; }

    /*!
     * \brief Set the DES Length Scale.
     */
    inline void SetDES_LengthScale(su2double val_des_lengthscale) {DES_LengthScale = val_des_lengthscale; }

    /*!
     * \brief Set the vortex tilting measure for computation of the EDDES length scale
     */
    void SetVortex_Tilting(su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity);

    /*!
     * \brief Get the vortex tilting measure for computation of the EDDES length scale
     * \return Value of the DES length Scale
     */
    inline su2double GetVortex_Tilting() {return Vortex_Tilting; }
    /*!
     * \brief Get the machine learning parameter.
     * \param[in] par_index - Index of point.
     * \return Value of the machine learning parameter.
     */
    su2double GetMLParam() override {return ML_Param; }
    /*!
     * \brief Set the machine learning parameter.
     * \param[in] par_index - Index of point.
     * \param[in] val_mlparam - New value of the machine learning parameter.
     */
    void SetMLParam(su2double val_mlparam) override {
        ML_Param = val_mlparam;
    }
    /*!
     * \brief Get the original index.
     * \param[in] par_index - Index of point.
     * \return Value of the machine learning parameter.
     */
    su2double GetMLIndex() override {return original_index; }
    /*!
     * \brief Set the original index.
     * \param[in] val_index - New value of the index.
     */
    void SetMLIndex(unsigned long val_index) override {
        original_index = val_index;
    }
};










