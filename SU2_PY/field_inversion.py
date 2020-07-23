#!/usr/bin/env python

## \file field_inversion.py
#  \Author R. Pochampalli
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
from shutil import copyfile

# -------------------------------------------------------------------
#  Auxiliary functions
# -------------------------------------------------------------------

def _isfinite(val):
    if math.isinf(val) or math.isnan(val):
        raise ValueError


# updates the parameters in the config files
def update_config(filename, params):
    fid = open(filename, "r")
    lines = fid.readlines()
    fid.close()

    for param in params:
        for i in range(len(lines)):
            if lines[i].startswith(param):
                lines[i] = param + "= " + repr(params[param]) + "\n"
                break
            # end
        # end
    # end

    fid = open(filename, "w")
    fid.writelines(lines)
    fid.close()
    # end




# end


# -------------------------------------------------------------------
#  Project class
# -------------------------------------------------------------------

class Project:
    def __init__(self, commands, inputFile, configFiles, outputFiles):
        self._inputFile = inputFile
        self._objValFile = outputFiles[0]
        self._objDerFile = outputFiles[1]
        self._objValCommand = commands[0] + configFiles[0] + " > objval.stdout"
        self._objDerCommand = commands[1] + configFiles[1] + " > objder.stdout"
        self._directConfigFile = configFiles[0]
        self._adjConfigFile = configFiles[1]
        self._sol_file = 'solution_flow.dat'
        self._res_file = 'restart_flow.dat'
        self._sol_file_adj = 'solution_adj.dat'
        self._res_file_adj = 'restart_adj.dat'
        self._val_reg_param = 1e-5
        self._niter = 0

    # end

    def obj_val(self, x, itn):
        obj_fun = 0.0
        self._niter = itn['Eval']
        if self._niter == 0:
            try:
                os.remove(self._objValFile)
            except:
                pass
            if Path('./' + self._sol_file).is_file():
                if Path('./' + self._res_file).is_file():
                    Path('./' + self._sol_file).unlink()
                    copyfile(self._res_file, self._sol_file)
                    print("Running a restart.")
                else:
                    print("Running cold: primal solver.")
            else:
                if Path('./' + self._res_file).is_file():
                    copyfile(self._res_file, self._sol_file)
                    print("Running a restart.")
                else:
                    print("Running cold: primal solver.")
        else:
            os.rename(self._objValFile, str(niter) + '_' + self._objValFile)
            os.rename(self._sol_file, str(niter) + '_' + self._sol_file)
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
        obj_fun = ofr[0] + self._val_reg_param * self.regularization(x)
        self._niter += 1

        return obj_fun

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
                    print("Running cold: adjoint solver.")
            if Path('./' + self._sol_file).is_file():
                if Path('./' + self._res_file).is_file():
                    Path('./' + self._sol_file).unlink()
                    os.rename(self._res_file, self._sol_file)
                else:
                    print("Solution file not found.")
            else:
                if Path('./' + self._res_file_adj).is_file():
                    os.rename(self._res_file_adj, self._sol_file_adj)
                else:
                    print("Running cold: adjoint solver.")
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
        grads = np.asarray(grads, dtype=np.double)
        grad_reg = self.grad_regularization(x)
        if grad_reg is not None:
            grads = np.add(grads, grad_reg)
        return grads

    #: obj_der()

    def update_params(self, param_upd):
        par_file = open(self._inputFile, 'a+')
        par_file.seek(0)
        par_file.truncate()
        par_file.close()

        with open(self._inputFile, 'a') as par_file:
            # par_file.write("NPARA=" + str(len(param_upd)) + "\n")
            for item in param_upd:
                par_file.write(str(item) + "\n")

    #: update_params()

    def regularization(self, params):
        val_reg = 0.0
        for val_param in params:
            val_reg += pow(val_param - 1.0, 2)
        val_reg = math.sqrt(val_reg)
        return val_reg

    def grad_regularization(self, params):
        val_reg = self.regularization(params)
        if val_reg == 0.0:
            return
        else:
            pd_obj_x = self._val_reg_param * (params - 1.0) / val_reg
            return pd_obj_x


#: Project


# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():
    # general options for L-BFGS-B
    options = {'disp': True, 'maxcor': 10, 'ftol': 1e-16, 'gtol': 1e-16, 'maxiter': 100, 'maxls': 20}

    # these are the commands for the direct and adjoint runs, modify to run parallel
    commands = ["mpirun SU2_CFD ", "mpirun SU2_CFD_AD "]

    # file through which SU2 gets the parameter values
    inputFile = "ml_param.su2"

    # names of the output files [objective value, objective gradient, constraint value, ...]
    outputFiles = ["of_func_file.dat", "field_sensitivity.csv"]

    # settings for direct run and adjoint of the objective
    fnames = ["direct.cfg", "adjoint.cfg"]

    params = np.loadtxt('ml_param.su2', dtype=np.double, delimiter='\n')
    # params = np.asarray([np.double(p_val) for p_val in open('inputFile', 'r') if p_val[0] != 'N'])
    obj = Project(commands, inputFile, fnames, outputFiles)
    line = "# -------------------------------------------------------------------\n"
    message = "#  Begin Optimization\n"
    print(line + message + line)
    opt = scipy.optimize.minimize(obj.obj_val, params, args=({'Eval': 0},), method="L-BFGS-B", jac=obj.obj_der,
                                  bounds=None, options=options)
    print(opt)


#: def main()

# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

if __name__ == '__main__':
    main()
