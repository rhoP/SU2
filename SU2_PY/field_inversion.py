#!/usr/bin/env python

## \file field_inversion.py
#  \brief Python script for field inversion.
#  \version 7.0.3 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#

import os
import math
import subprocess as sp
import numpy as np
import scipy.optimize
from pathlib import Path

# -------------------------------------------------------------------
#  Project class
# -------------------------------------------------------------------

def _isfinite(val):
    if math.isinf(val) or math.isnan(val):
        raise ValueError


# updates the parameters in the config files
def update_config(fnames, params):
    for fname in fnames:
        fid = open(fname, "r")
        lines = fid.readlines()
        fid.close()

        for param in params:
            for i in range(len(lines)):
                if lines[i].startswith(param.name()):
                    lines[i] = param.name() + "= " + repr(param.value()) + "\n"
                    break
                # end
            # end
        # end

        fid = open(fname, "w")
        fid.writelines(lines)
        fid.close()
    # end


# end


class Project:
    def __init__(self, commands, inputFile, configFiles, outputFiles):
        self._inputFile = inputFile
        self._objValFile = outputFiles[0]
        self._objDerFile = outputFiles[1]
        self._objValCommand = commands[0] + configFiles[0] + " > objval.stdout"
        self._objDerCommand = commands[1] + configFiles[1] + " > objder.stdout"
        self._sol_file = 'solution_flow.dat'
        self._res_file = 'restart_flow.dat'
        self._sol_file_adj = 'solution_adj_turbML.dat'
        self._res_file_adj = 'restart_adj_turbML.dat'
    # end

    def obj_val(self, x, itn):
        if itn['Eval'] % 100 == 0:
            try:
                os.remove(self._objValFile)
            except:
                pass
            if Path('./'+self._sol_file).is_file():
                if Path('./'+self._res_file).is_file():
                    Path('./' + self._sol_file).unlink()
                    os.rename(self._res_file, self._sol_file)
                else:
                    print("Running cold.")
            else:
                if Path('./'+self._res_file).is_file():
                    os.rename(self._res_file, self._sol_file)
                else:
                    print("Running cold.")
        else:
            niter = itn['Eval']
            os.rename(self._objValFile, str(niter)+'_'+self._objValFile)
            os.rename(self._sol_file, str(niter)+'_'+self._sol_file)
        try:
            self.update_params(x)
            sp.call(self._objValCommand, shell=True)
            ofr = [np.double(g_val) for g_val in open(self._objValFile, 'r')]
            if len(ofr) > 1:
                raise RuntimeError("Objective function evaluation failed: too many inputs.")
            _isfinite(ofr[0])
        except:
            raise RuntimeError("Objective function evaluation failed")
        # end
        return ofr[0]

    #: Obj_val()

    def obj_der(self, x, itn):
        if itn['Eval'] % 100 == 0:
            try:
                os.remove(self._objDerFile)
            except:
                pass
            if Path('./' + self._sol_file_adj).is_file():
                if Path('./' + self._res_file_adj).is_file():
                    Path('./' + self._sol_file_adj).unlink()
                    os.rename(self._res_file_adj, self._sol_file_adj)
                else:
                    print("Running cold.")
            else:
                if Path('./' + self._res_file_adj).is_file():
                    os.rename(self._res_file_adj, self._sol_file_adj)
                else:
                    print("Running cold.")
        else:
            niter = itn['Eval']
            os.rename(self._objDerFile, str(niter) + '_' + self._objDerFile)
            os.rename(self._sol_file_adj, str(niter) + '_' + self._sol_file_adj)
        try:
            # main command
            sp.call(self._objDerCommand, shell=True)
            grads = [np.double(g_val) for g_val in open(self._objDerFile)]
            if len(grads) != x.shape[0]:
                raise RuntimeError("Mismatch in dimension of gradients")

            for val in grads:
                _isfinite(val)
            # end
        except:
            raise RuntimeError("Gradient evaluation failed")
        # end
        return np.asarray(grads, dtype=np.double)

    #: obj_der()

    def update_params(self, param_upd):
        par_file = open(self._inputFile, 'a+')
        par_file.seek(0)
        par_file.truncate()
        par_file.close()

        with open(self._inputFile, 'a') as par_file:
            par_file.write("NPARA=" + str(len(param_upd)) + "\n")
            for item in param_upd:
                par_file.write(str(item) + "\n")

    #: update_params()


#: Project


# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():
    # general options for L-BFGS-B
    options = {'disp': True, 'maxcor': 10, 'ftol': 1e-5, 'gtol': 1e-18, 'maxiter': 100, 'maxls': 2}

    # these are the commands for the direct and adjoint runs, modify to run parallel
    commands = ["mpirun -np 6 --oversubscribe SU2_CFD ", "SU2_CFD_AD "]

    # file through which SU2 gets the parameter values
    inputFile = "ml_param.su2"

    # names of the output files [objective value, objective gradient, constraint value, ...]
    outputFiles = ["of_func_file.dat", "field_sensitivity.csv"]

    # settings for direct run and adjoint of the objective
    fnames = ["direct.cfg", "adjoint.cfg"]

    params = np.full(shape=14576, fill_value=1.0, dtype=np.double)
    # params = np.asarray([np.double(p_val) for p_val in open('inputFile', 'r') if p_val[0] != 'N'])
    obj = Project(commands, inputFile, fnames, outputFiles)
    line = "# -------------------------------------------------------------------\n"
    message = "#  Begin Optimization\n"
    print(line + message + line)
    opt = scipy.optimize.minimize(obj.obj_val, params, args=({'Eval':0},),method="L-BFGS-B", jac=obj.obj_der,
                                  bounds=None, options=options)


#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
