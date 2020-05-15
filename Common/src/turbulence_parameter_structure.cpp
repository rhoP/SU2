
/*!
 * \file CMeshReaderFVM.hpp
 * \brief Header file for the class CMeshReaderFVM.
 *        The implementations are in the <i>CMeshReaderFVM.cpp</i> file.
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

#include "../include/turbulence_parameter_structure.hpp"
#include "../include/toolboxes/CLinearPartitioner.hpp"

#include <iostream>
#include <cstdlib>

/*
     * \brief  Constructor of the class CTurbML
     * \param[in] nParam - Number of Parameters of the Turbulence model.
     * \param[in] config - Definition of the particular problem.
*/
CTurbML::CTurbML(CConfig  *val_config, unsigned short val_iZone,
                 unsigned short val_nZone)
        :CMeshReaderFVM(val_config, val_iZone, val_nZone){

    this->config = val_config;


    /* Store the current zone to be read and the total number of zones. */
    myZone = val_iZone;
    nZones = val_nZone;


    /* Store the mesh filename since we will open/close multiple times. */

    MLParam_Filename = config->GetMLParam_FileName();
    /* Read the basic metadata and perform some basic error checks. */

    ReadMetadata();



    /* Read and store the parameter values */

    ReadParameterValues();

}



void CTurbML::ReadMetadata() {
        MLParam_file.open(MLParam_Filename.c_str(), ios::in);

                /*--- Check if parameter file is open ---*/
                        if (MLParam_file.fail()) {
                SU2_MPI::Error(string("Error opening parameter file.") +
                                                       string(" \n Check if the file exists."), CURRENT_FUNCTION);
            }
        string text_line;
        string::size_type position;
        /*--- Read the metadata: total number of machine learning parameters. ---*/

        bool foundNPARA = false;

        while (getline (MLParam_file, text_line)) {

            /*--- Read the number of parameters of the problem ---*/

            position = text_line.find ("NPARA=",0);
            if (position != string::npos) {
                text_line.erase (0,6);
                numberOfMLParameters = atoi(text_line.c_str());
                for (unsigned long iPara = 0; iPara < numberOfMLParameters; iPara++)
                    getline (MLParam_file, text_line);
                foundNPARA = true;
            }
        }

                /* Throw an error if the parameter keyword was not found. */
                        if (!foundNPARA) {
                SU2_MPI::Error(string("Could not find NPARA= keyword.") +
                                                       string(" \n Check the SU2 parameter file format."),
                                +                       CURRENT_FUNCTION);
            }
        MLParam_file.close();
    }

void CTurbML::ReadParameterValues() {


    /*--- Reserve memory for the vector of parameters ---*/
    ML_Parameters.reserve(numberOfMLParameters);
    MLParam_file.open(MLParam_Filename.c_str(), ios::in);

    /*--- Read the parameters into our data structure. ---*/
    string text_line;
    string::size_type position;


    while (getline (MLParam_file, text_line)) {

        position = text_line.find("NPARA=",0);
        if (position != string::npos) {

            unsigned long GlobalIndex = 0;
            while (GlobalIndex < numberOfMLParameters) {

                getline(MLParam_file, text_line);

                /*--- We only read information for this node if it is owned by this
                 rank based upon our initial linear partitioning. ---*/

                    double par_val{0.0};
                    istringstream par_line(text_line);
                    par_line >> par_val;
                    ML_Parameters.emplace_back(par_val);

                GlobalIndex++;
            }
        }
    }

    MLParam_file.close();
}


bool CTurbML::MatchParamsPoints(unsigned long global_points) const {
    return numberOfMLParameters == global_points;
}