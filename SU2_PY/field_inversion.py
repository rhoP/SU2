#!/usr/bin/env python

## \file field_inversion.py
#  \brief Python script for performing the shape optimization.
#  \author R. Pochampalli
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

import os, sys, shutil
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():

    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-r", "--name", dest="projectname", default='',
                      help="try to restart from project file NAME", metavar="NAME")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="DISCRETE_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-o", "--optimization", dest="optimization", default="SLSQP",
                      help="OPTIMIZATION techique (SLSQP, CG, BFGS, POWELL)", metavar="OPTIMIZATION")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")


    (options, args)=parser.parse_args()

    # process inputs
    options.partitions  = int( options.partitions )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    options.nzones      = int( options.nzones )

    sys.stdout.write('\n-----------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )   Release 7.0.3 \"Blackbird\"                       |\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                    |\n')
    sys.stdout.write('|   |___/\\___//___|  Field Inversion for Turbulence Modeling           |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Project Website: https://su2code.github.io                        |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| The SU2 Project is maintained by the SU2 Foundation                   |\n')
    sys.stdout.write('| (http://su2foundation.org)                                            |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)                |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is free software; you can redistribute it and/or                  |\n')
    sys.stdout.write('| modify it under the terms of the GNU Lesser General Public            |\n')
    sys.stdout.write('| License as published by the Free Software Foundation; either          |\n')
    sys.stdout.write('| version 2.1 of the License, or (at your option) any later version.    |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is distributed in the hope that it will be useful,                |\n')
    sys.stdout.write('| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n')
    sys.stdout.write('| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n')
    sys.stdout.write('| Lesser General Public License for more details.                       |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| You should have received a copy of the GNU Lesser General Public      |\n')
    sys.stdout.write('| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')

    field_inversion( options.filename    ,
                        options.projectname ,
                        options.partitions  ,
                        options.gradient    ,
                        options.optimization ,
                        options.quiet       ,
                        options.nzones      )

#: main()

def field_inversion   ( filename                           ,
                        projectname  = ''                  ,
                        partitions   = 0                   ,
                        gradient     = 'DISCRETE_ADJOINT'  ,
                        optimization = 'SLSQP'             ,
                        quiet        = False               ,
                        nzones       = 1                   ):
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )
    if quiet: config.CONSOLE = 'CONCISE'
    config.GRADIENT_METHOD = gradient

    its               = int ( config.OPT_ITERATIONS )                      # number of opt iterations
    bound_upper       = float ( config.OPT_BOUND_UPPER )                   # variable bound to be scaled by the line search
    bound_lower       = float ( config.OPT_BOUND_LOWER )                   # variable bound to be scaled by the line search
    relax_factor      = float ( config.OPT_RELAX_FACTOR )                  # line search scale
    gradient_factor   = float ( config.OPT_GRADIENT_FACTOR )               # objective function and gradient scale
    def_dv            = config.DEFINITION_DV                               # complete definition of the desing variable #TODO review
    n_dv              = sum(def_dv['SIZE'])                                # number of design variables                 #TODO review
    accu              = float ( config.OPT_ACCURACY ) * gradient_factor    # optimizer accuracy
    x0                = [0.0]*n_dv                                         # initial design                             #TODO review
    xb_low            = [float(bound_lower)/float(relax_factor)]*n_dv      # lower dv bound it includes the line search acceleration factor
    xb_up             = [float(bound_upper)/float(relax_factor)]*n_dv      # upper dv bound it includes the line search acceleration fa
    xb                = list(zip(xb_low, xb_up)) # design bounds

    # State
    state = SU2.io.State()
    state.find_files(config)

    # Project
    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project.config = config
    else:
        project = SU2.opt.Project(config,state)

    # Optimize
    if optimization == 'SLSQP':
      SU2.opt.SLSQP(project,x0,xb,its,accu)
    if optimization == 'CG':
      SU2.opt.CG(project,x0,xb,its,accu)
    if optimization == 'BFGS':
      SU2.opt.BFGS(project,x0,xb,its,accu)
    if optimization == 'POWELL':
      SU2.opt.POWELL(project,x0,xb,its,accu)


    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)

    return project

#: field_inversion()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
