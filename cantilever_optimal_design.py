###############################################################################################
#
# Author: Jack Adams
# Date: 9/01/2018
#
# Optimal Design Script
#
# This script is used to predict the optimal combination of design variables for an experiment.
# This specifically will work for a cantilever beam fixed at one end.
#
# The design variables in question are two angles, theta and phi, which represent the amount
# the beam is lifted from the horizontal and the angle which it is twisted around its
# longitudinal axis by.
#
# This script assumes that the input parameter value is the true material parameter, even if
# there is no evidence for this being the case. This method will assume that the value is true
# and then simulate the deformation at a series of angle combinations. At each, the value of
# the Hessian matrix will be found.
#
# This matrix will provide information about the determinability of the material parameter,
# namely how sensitive that configuration is to changes in the material parameter. The more
# sensitive the angles combination is, the more likely the results will provide an accurate
# measure of the material parameter.
#
##############################################################################################

from parameter_optimisation import ParameterEstimation
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
import numpy as np
import math

# First set up the variables needed to create the simulation.
dimensions = np.array([30, 12, 12])
elements = np.array([2, 1, 2])
parameter_value = np.array([7.0, 1.5])

# Next, create the instance of the simulation class and add the initialised variables to it.
ps = ParameterEstimation()
ps.simulation = CantileverSimulation()
ps.simulation.set_cantilever_dimensions(dimensions)
ps.simulation.set_cantilever_elements(elements)
#ps.simulation.set_Xi_points_num(4)
ps.simulation.set_diagnostic_level(1)
ps.simulation.setup_cantilever_simulation()
ps.simulation.set_Neo_Hookean_single_layer(parameter_value)

# Now define the design variables.
thetaStart = -90
thetaEnd = 90
thetaStep = 30
phiStart = 0
phiEnd = 180
phiStep = 30

detHMatrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))
condHMatrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))
detH0Matrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))

loopCounter = 1
loopMax = (((thetaEnd - thetaStart) / thetaStep) + 1) * (((phiEnd - phiStart) / phiStep) + 1)

# Now loop through the design variables and solve the simulation under each condition.
for theta in range(thetaStart, thetaEnd+1, thetaStep):
    for phi in range(phiStart, phiEnd+1, phiStep):
        grav_vect = ps.simulation.gravity_vector_calculation((theta * math.pi / 180), (phi * math.pi / 180))
        ps.simulation.set_gravity_vector(grav_vect)
        ps.simulation.solve_simulation()
        ps.simulation.set_projection_data()

        # Next calculate the Hessian matrix for each design variable combination.
        ps.set_objective_function(two_layer_objective_function)
        [H, detH, condH, detH0] = ps.evaluate_hessian(parameter_value, 1e-7)

        print "Simulation {0} of {1}: Complete.".format(loopCounter, loopMax)
        print 'Determinant of Hessian = {0}'.format(detH)
        print '\n'

        loopCounter += 1

        # Now compile these into a matrix
        detHMatrix[theta/thetaStep + 1, phi/phiStep] = detH
        condHMatrix[theta/thetaStep + 1, phi/phiStep] = condH
        detH0Matrix[theta/thetaStep + 1, phi/phiStep] = detH0


print 'Optimal Design finished.'
print '\n'
print 'Printing results to file for visualisation:'

# Export the matrix results for visualisation.
np.savetxt('optimal_design_detH.txt', detHMatrix, delimiter=' ', newline='\n')
designVariables = np.array([thetaStart, thetaEnd, thetaStep, phiStart, phiEnd, phiStep])
np.savetxt('optimal_design_variables.txt', designVariables, delimiter=" ", newline="\n")

print 'Files generated'
