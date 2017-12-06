import numpy as np
import CantileverSimulation as CantileverSimulation
from CantileverSimulation import cantilever_objective_function
import ParameterOptimisation as ParamOpt
import math

# First gather all the data sets. If this is artificial, run a simulation with the known parameters first and extract
# points to use first. If this is real data, read in all the data before starting the loops and just access the correct
# orientation's data in each loop.

# Start by setting some parameters and creating the simulation.
cantilever_dimensions = np.array([60, 40, 40])
cantilever_elements = np.array([1, 1, 1])
true_parameter = np.array([1.452])
guess_parameter = np.array([0.5])

ps = ParamOpt()
ps.simulation = CantileverSimulation()

# Now enter into the loop and set up each simulation with a different gravity vector.
designVariableOneStart = -90
designVariableOneFinish = 90
designVariableOneStep = 15
designVariableTwoStart = 0
designVariableTwoFinish = 360
designVariableTwoStep = 20

HMatrix = detHMatrix = condHMatrix = detH0Matrix = np.zeros(((designVariableOneFinish-designVariableOneStart)/designVariableOneStep, (designVariableTwoFinish-designVariableTwoStart)/designVariableTwoStep))

for i in range(designVariableOneStart, designVariableOneFinish, designVariableOneStep):
    for j in range(designVariableTwoStart, designVariableTwoFinish, designVariableTwoStep):
        gravity_vector = ps.simulation.gravity_vector_calculation(i*math.pi/180, j*math.pi/180)

        ps.simulation.set_cantilever_dimensions(cantilever_dimensions)
        ps.simulation.set_cantilever_elements(cantilever_elements)
        ps.simulation.set_gravity_vector(gravity_vector)
        ps.simulation.set_diagnostic_level(0)
        ps.simulation.setup_cantilever_simulation()
        ps.simulation.prepare_projection()
        ps.simulation.set_Mooney_Rivlin_parameter_values(true_parameter)
        ps.simulation.solve_simulation()
        data = ps.simulation.generate_data(3)

        simulation_tuple = (ps.simulation,)
        ps.set_objective_function(cantilever_objective_function, simulation_tuple)

# Now enter into the loops which change the design variables and

# Next, run FE models with a guess as to what the material property is and optimise the parameter to provide the best
# deformation match to the provided data.


# Now that the optimal parameters have been found for this orientation, calculate the Hessian and its properties in this
# orientation and store the results.


# Once all the orientations have been analysed, plot the resulting Hessian metrics against the orientation design
# variables.
