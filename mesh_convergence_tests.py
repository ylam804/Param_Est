####################################################################################################################
#
# Author: Jack Adams
# Date: 2017/12/05
#
#
# mesh_convergence_tests:
#
# This script will run a finite element model and data creation routine to determine how many elements are required
# to reach mesh independence, that is, where the elements in the mesh have no influence on the results of the
# simulation.
#
# This is done by repeatedly calling a data creation function after the FE model has been run which will determine
# the location of certain points in the model. It will then check if these location values have converged to within
# a tolerance, showing that the mesh no longer has an effect on the results of the simulation.
#
# The output will be a number of elements which are required to reach this point, specifying the number of elements
# in each direction.
#
####################################################################################################################

import numpy as np
import CantileverSimulation
import ParameterOptimisation



