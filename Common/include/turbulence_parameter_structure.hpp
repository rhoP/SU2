
/*!
 * \file turbulence_parameter_structure.hpp
 * \brief Header file for the class CTurbML.
 *        The implementations are in the <i>turbulence_parameter_structure.cpp</i> file.
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

#include "./mpi_structure.hpp"
#include "./CConfig.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>


/*!
 * \class CTurbML
 * \brief Main class for managing the parameters for machine learning with turbulence models.
 * \ingroup Turbulence_Model
 * \author Rohit P.
 */

class CTurbML {

private:
    CConfig *config = nullptr;                         /*!< \brief Local pointer to the config parameter object. */
    unsigned long numberOfMLParameters;      /*!< \brief Number of parameter values in the parameter file. */
    string MLParam_Filename;                 /*!< \brief Name of the SU2 Parameter file being read. */
    ifstream MLParam_file;                   /*!< \brief File object for the SU2 ASCII mesh file. */
    std::vector<su2double> ML_Parameters;    /*!< \brief Vector containing the parameter values. */
    /*!
     * \brief Reads all SU2 ASCII mesh metadata and checks for errors.
     */
    void ReadMetadata();
    /*!
     * \brief Reads the grid points from an SU2 zone into linear partitions across all ranks.
     */
    void ReadParameterValues();
public:
    /*!
     * \brief Constructor of the CMLParamReader class.
     */
    CTurbML(CConfig *val_config, unsigned long global_points);
    /*!
     * \brief Destructor of the CMLParamReader class.
     */
    ~CTurbML() = default;
    /*!
     * \brief Get the machine learning parameter.
     * \param[in] par_index - Index of point.
     * \return Value of the machine learning parameter.
     */
    su2double Get_iParamML(unsigned long par_index) {return ML_Parameters[par_index]; }
    /*!
     * \brief Set the machine learning parameter.
     * \param[in] par_index - Index of point.
     * \param[in] val_mlparam - New value of the machine learning parameter.
     */
    void Set_iParamML(su2double val_mlparam, unsigned long par_index) {
        ML_Parameters[par_index] = val_mlparam;
    }
    /*!
     * \brief Get the number of machine learning parameters.
     * \return Value of the machine learning parameter.
     */
    unsigned long Get_nParamML() {return numberOfMLParameters; }

    /*!
     * \brief Set the number of machine learning parameters.
     * \return Value of the machine learning parameter.
     */
    void Set_nParamsML(unsigned long val_nParams) {numberOfMLParameters = val_nParams; }

    /*!
     * \brief Match the number of parameters with the global number of points.
     * \return global number of points.
     */
    void MatchParamsPoints(unsigned long global_points);
/*!
     * \brief Get the ith index from the RCM ordering vector.
     * \return index of original ordering.
     */
    std::vector<unsigned long> RCM_ordering; /*!< \brief Vector containing result vector of the RCM ordering. */
    unsigned long Get_iRCM_Result(unsigned long ind) {return RCM_ordering[ind]; }

    /*!
     * \brief Set the number of machine learning parameters.
     * \return Value of the machine learning parameter.
     */
    void Set_iRCM_Order(unsigned long ind, unsigned long val_ind) {RCM_ordering[ind] = val_ind; }

};




























