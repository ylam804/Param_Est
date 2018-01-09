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
from cantilever_simulation import cantilever_objective_function
import numpy as np
import math

# First set up the variables needed to create the simulation.
dimensions = np.array([30, 12, 12])
elements = np.array([2, 1, 1])
parameter_value = np.array([2.05])

# Next, create the instance of the simulation class and add the initialised variables to it.
ps = ParameterEstimation()
ps.simulation = CantileverSimulation()
ps.simulation.set_cantilever_dimensions(dimensions)
ps.simulation.set_cantilever_elements(elements)
ps.simulation.set_Xi_points_num(4)
ps.simulation.set_diagnostic_level(0)
ps.simulation.setup_cantilever_simulation()
ps.simulation.set_Mooney_Rivlin_parameter_values(parameter_value)

# Now define the design variables.
#thetaStart = -90
#thetaEnd = 90
#thetaStep = 90
#phiStart = 0
#phiEnd = 180
#phiStep = 90

theta = 45
phi = 20

#HMatrix = detHMatrix = condHMatrix = detH0Matrix = np.zeros((((thetaEnd - thetaStart) / thetaStep) + 1, ((phiEnd - phiStart) / phiStep) + 1))

# Now loop through the design variables and solve the simulation under each condition.
#for theta in range(thetaStart, thetaEnd+1, thetaStep):
#    for phi in range(phiStart, phiEnd+1, phiStep):
grav_vect = ps.simulation.gravity_vector_calculation((theta * math.pi / 180), (phi * math.pi / 180))
ps.simulation.set_gravity_vector(grav_vect)
ps.simulation.solve_simulation()
ps.simulation.set_projection_data()

ps.set_objective_function(cantilever_objective_function)
[H, detH, condH, detH0] = ps.evaluate_hessian(parameter_value, 1e-7)

print detH

#        HMatrix[theta/thetaStep, phi/phiStep] = H
#        detHMatrix[theta/thetaStep, phi/phiStep] = detH
#        condHMatrix[theta/thetaStep, phi/phiStep] = condH
#        detH0Matrix[theta/thetaStep, phi/phiStep] = detH0

# Next calculate the Hessian matrix for each design variable combination.
#print detHMatrix

# Now compile these into a matrix and export them for visualisation.
