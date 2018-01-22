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

def destroy_routine(simulation):
    """
    Destroys some parts of the simulation so that they can be re-initialised with new values.
    """

    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()

dimensions = np.array([30, 12, 12])
elements = np.array([3, 3, 3])
initial_parameter = np.array([8.4378])

# Now define the design variables.
thetaStart = -90
thetaEnd = 90
thetaStep = 45
phiStart = 0
phiEnd = -179
phiStep = -45

theta = np.array((range(thetaStart, thetaEnd+1, thetaStep))) * math.pi / 180
phi = np.array((range(phiStart, phiEnd+1, phiStep))) * math.pi / 180
detHMatrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))
condHMatrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))
detH0Matrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))

loopCounter = 1
loopMax = (((thetaEnd - thetaStart) / thetaStep) + 1) * (((phiEnd - phiStart) / phiStep) + 1)

# Now loop through the design variables and solve the simulation under each condition.
for i in range(len(theta)):
    for j in range(len(phi)):

        ps = ParameterEstimation()
        ps.simulation = CantileverSimulation()
        ps.simulation.set_cantilever_dimensions(dimensions)
        ps.simulation.set_cantilever_elements(elements)
        ps.simulation.set_gravity_vector(ps.simulation.gravity_vector_calculation(theta[i], phi[j]))
        ps.simulation.set_diagnostic_level(0)
        ps.simulation.setup_cantilever_simulation()
        ps.simulation.set_Neo_Hookean_single_layer(initial_parameter)
        ps.simulation.solve_simulation()

        #ps.simulation.export_results("Cantilever")
        ps.simulation.set_projection_data()
        ps.set_objective_function(single_layer_objective_function)
        [H, detH, condH, detH0] = ps.new_evaluate_hessian_method(initial_parameter, 1e-7)

        print "Simulation {0} of {1}: Complete.".format(loopCounter, loopMax)
        print "For angles Theta = {0}, Phi = {1}".format(theta[i], phi[j])
        print "     Gravity X-Component = {0}".format(ps.simulation.gravity_vector[0])
        print "     Gravity Y-Component = {0}".format(ps.simulation.gravity_vector[1])
        print "     Gravity Z-Component = {0}".format(ps.simulation.gravity_vector[2])
        print "Determinant of Hessian = {0}".format(detH)
        print "\n"

        loopCounter += 1

        # Now compile these into a matrix
        detHMatrix[i, j] = detH
        condHMatrix[i, j] = condH
        detH0Matrix[i, j] = detH0

        #destroy_routine(ps.simulation)
        ps = None


print 'Optimal Design finished.'
print '\n'
print 'Printing results to file for visualisation:'

# Export the matrix results for visualisation.
np.savetxt('optimal_design_detH.txt', detHMatrix, delimiter=' ', newline='\n')
designVariables = np.array([thetaStart, thetaEnd, thetaStep, phiStart, phiEnd, phiStep])
np.savetxt('optimal_design_variables.txt', designVariables, delimiter=" ", newline="\n")

print 'Files generated'
