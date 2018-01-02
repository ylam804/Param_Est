import numpy as np
from opencmiss.iron._utils import CMFEError
from CantileverSimulation import CantileverSimulation
from CantileverSimulation import cantilever_objective_function
from CantileverSimulation import destroy_routine
from ParameterOptimisation import ParameterEstimation
import math

# First gather all the data sets. If this is artificial, run a simulation with the known parameters first and extract
# points to use first. If this is real data, read in all the data before starting the loops and just access the correct
# orientation's data in each loop.

# Start by setting some parameters and creating the simulation.
cantilever_dimensions = np.array([60, 40, 40])
cantilever_elements = np.array([2, 2, 2])
true_parameter = np.array([1.452])
guess_parameter = np.array([0.5])

ps = ParameterEstimation()
ps.simulation = CantileverSimulation()

# Now enter into the loop and set up each simulation with a different gravity vector.
designVariableOneStart = -90
designVariableOneFinish = 90
designVariableOneStep = 30
designVariableTwoStart = 0
designVariableTwoFinish = 360
designVariableTwoStep = 45

angleOneCounter = 0
angleTwoCounter = 0

HMatrix = detHMatrix = condHMatrix = detH0Matrix = np.zeros((((designVariableOneFinish-designVariableOneStart)/designVariableOneStep) + 1, ((designVariableTwoFinish-designVariableTwoStart)/designVariableTwoStep) + 1))

for i in range(designVariableOneStart, designVariableOneFinish+1, designVariableOneStep):
    for j in range(designVariableTwoStart, designVariableTwoFinish+1, designVariableTwoStep):
        gravity_vector = ps.simulation.gravity_vector_calculation(i*math.pi/180, j*math.pi/180)

        ps.simulation.set_cantilever_dimensions(cantilever_dimensions)
        ps.simulation.set_cantilever_elements(cantilever_elements)
        ps.simulation.set_gravity_vector(gravity_vector)
        ps.simulation.set_diagnostic_level(0)
        ps.simulation.setup_cantilever_simulation()
        ps.simulation.set_Mooney_Rivlin_parameter_values(true_parameter)

        try:
            ps.simulation.solve_simulation()
        except CMFEError:
            [HMatrix[angleOneCounter,angleTwoCounter], detHMatrix[angleOneCounter,angleTwoCounter], condHMatrix[angleOneCounter,angleTwoCounter], detH0Matrix[angleOneCounter,angleTwoCounter]] = float('nan')

        if np.isnan(HMatrix[angleOneCounter,angleTwoCounter]) == False:

            data = ps.simulation.generate_data(3)
            ps.simulation.set_projection_data(data)
            ps.simulation.prepare_projection()

            # Now that the data points have been obtained, use them to find a set of optimised parameters.
            destroy_routine(ps.simulation)
            ps.initial_parameters = guess_parameter
            ps.simulation.set_projection_data(data)
            ps.simulation.setup_cantilever_simulation()
            ps.simulation.prepare_projection()
            simulation_tuple = (ps.simulation,)
            ps.set_objective_function(cantilever_objective_function, simulation_tuple)
            try:
                ps.optimise()
            except CMFEError:
                [HMatrix[angleOneCounter,angleTwoCounter], detHMatrix[angleOneCounter,angleTwoCounter], condHMatrix[angleOneCounter,angleTwoCounter], detH0Matrix[angleOneCounter,angleTwoCounter]] = [float('nan'), float('nan'), float('nan'), float('nan')]
            else:
                # Now use those optimised parameters to calculate the Hessian metrics.
                [HMatrix[angleOneCounter,angleTwoCounter], detHMatrix[angleOneCounter,angleTwoCounter], condHMatrix[angleOneCounter,angleTwoCounter], detH0Matrix[angleOneCounter,angleTwoCounter]] = ps.evaluate_hessian(ps.solutions.x, 1e-7)

        destroy_routine(ps.simulation)
        angleTwoCounter += 1
    angleOneCounter += 1
    angleTwoCounter = 0

# Once all the orientations have been analysed, plot the resulting Hessian metrics against the orientation design
# variables.
print 'Done!'
