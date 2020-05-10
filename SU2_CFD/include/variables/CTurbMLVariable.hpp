/*!
 * \file CTurbMLVariable.hpp
 * \brief Declaration of the variables of the SA turbulence model with machine learning.
 * \author F. Palacios, T. Economon
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

#include "CTurbVariable.hpp"

/*!
 * \class CTurbMLVariable
 * \brief Main class for defining the variables of turbulence modeling with machine learning.
 * \ingroup Turbulence_Model
 * \author R. Pochampalli
 */

class CTurbMLVariable final : public CTurbVariable {

private:
    VectorType field_param;         /*!< \brief Value of the field parameter for turbulence modeling. */
    VectorType gamma_BC;
    VectorType DES_LengthScale;
    VectorType Vortex_Tilting;

public:
    /*!
     * \brief Constructor of the class.
     * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
     * \param[in] val_muT  - The eddy viscosity
     * \param[in] npoint - Number of points/nodes/vertices in the domain.
     * \param[in] ndim - Number of dimensions of the problem.
     * \param[in] nvar - Number of variables of the problem.
     * \param[in] constants -
     * \param[in] config - Definition of the particular problem.
     */
    CTurbMLVariable(su2double val_nu_tilde, su2double val_muT, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

    /*!
     * \brief Destructor of the class.
     */
    ~CTurbMLVariable() override = default;

    /*!
     * \brief Set the harmonic balance source term.
     * \param[in] iPoint - Point index.
     * \param[in] iVar - Index of the variable.
     * \param[in] source - Value of the harmonic balance source term. for the index <i>iVar</i>.
     */
    inline void SetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar, su2double source) override { HB_Source(iPoint,iVar) = source; }

    /*!
     * \brief Get the harmonic balance source term.
     * \param[in] iPoint - Point index.
     * \param[in] iVar - Index of the variable.
     * \return Value of the harmonic balance source term for the index <i>val_var</i>.
     */
    inline su2double GetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar) const override { return HB_Source(iPoint,iVar); }

    /*!
     * \brief Get the intermittency of the BC transition model.
     * \param[in] iPoint - Point index.
     * \return Value of the intermittency of the BC transition model.
     */
    inline su2double GetGammaBC(unsigned long iPoint) const override { return gamma_BC(iPoint); }

    /*!
     * \brief Set the intermittency of the BC transition model.
     * \param[in] iPoint - Point index.
     * \param[in] val_gamma - New value of the intermittency.
     */
    inline void SetGammaBC(unsigned long iPoint, su2double val_gamma) override { gamma_BC(iPoint) = val_gamma; }

    /*!
     * \brief Get the DES length scale
     * \param[in] iPoint - Point index.
     * \return Value of the DES length Scale.
     */
    inline su2double GetDES_LengthScale(unsigned long iPoint) const override { return DES_LengthScale(iPoint); }

    /*!
     * \brief Set the DES Length Scale.
     * \param[in] iPoint - Point index.
     */
    inline void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) override { DES_LengthScale(iPoint) = val_des_lengthscale; }

    /*!
     * \brief Set the vortex tilting measure for computation of the EDDES length scale
     * \param[in] iPoint - Point index.
     */
    void SetVortex_Tilting(unsigned long iPoint, const su2double* const* PrimGrad_Flow,
                           const su2double* Vorticity, su2double LaminarViscosity) override;

    /*!
     * \brief Get the vortex tilting measure for computation of the EDDES length scale
     * \param[in] iPoint - Point index.
     * \return Value of the DES length Scale
     */
    inline su2double GetVortex_Tilting(unsigned long iPoint) const override { return Vortex_Tilting(iPoint); }


    /*!
   * \brief Get the field parameter for turbulence modeling
   * \return value of the parameter
   */
    su2double Get_FieldParam(unsigned long iPoint) override {
        return field_param(iPoint);}
    /*!
     * \brief Set the field parameter.
     * \param[in] par_index - Index of point.
     * \param[in] val_param - New value of the field parameter.
     */
    void Set_FieldParam(unsigned long iPoint, su2double val_param) override {
        field_param(iPoint) = val_param;
    }
};
