import numpy as np
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import cantilever_objective_function
from parameter_optimisation import ParameterEstimation
import math

# First gather all the data sets. If this is artificial, run a simulation with the known parameters first and extract
# points to use first. If this is real data, read in all the data before starting the loops and just access the correct
# orientation's data in each loop.

# Start by setting some parameters and creating the simulation.
cantilever_dimensions = np.array([30, 12, 12])
cantilever_elements = np.array([1, 1, 1])
true_parameter = np.array([1.452])

ps = ParameterEstimation()
ps.simulation = CantileverSimulation()
ps.simulation.set_cantilever_dimensions(cantilever_dimensions)
ps.simulation.set_cantilever_elements(cantilever_elements)
ps.simulation.set_diagnostic_level(0)
ps.simulation.setup_cantilever_simulation()
ps.simulation.set_Mooney_Rivlin_parameter_values(true_parameter)

# Now enter into the loop and set up each simulation with a different gravity vector.
designVariableOneStart = -90
designVariableOneFinish = 90
designVariableOneStep = 60
designVariableTwoStart = 0
designVariableTwoFinish = 90
designVariableTwoStep = 45

angleOneCounter = 0
angleTwoCounter = 0

detHMatrix = np.zeros((((designVariableOneFinish-designVariableOneStart)/designVariableOneStep) + 1, ((designVariableTwoFinish-designVariableTwoStart)/designVariableTwoStep) + 1))
HMatrix = condHMatrix = detH0Matrix = np.zeros((((designVariableOneFinish-designVariableOneStart)/designVariableOneStep) + 1, ((designVariableTwoFinish-designVariableTwoStart)/designVariableTwoStep) + 1))

for i in range(designVariableOneStart, designVariableOneFinish+1, designVariableOneStep):

    print('\n')
    print angleOneCounter
    print('\n\n\n')

    for j in range(designVariableTwoStart, designVariableTwoFinish+1, designVariableTwoStep):

        gravity_vector = ps.simulation.gravity_vector_calculation(i*math.pi/180, j*math.pi/180)
        ps.simulation.set_gravity_vector(gravity_vector)
        ps.simulation.solve_simulation()
        data = ps.simulation.generate_data(3)
        ps.simulation.set_projection_data(data)
        simulation_tuple = (ps.simulation, )

        ps.set_objective_function(cantilever_objective_function, simulation_tuple)
        [H, detH, condH, detH0] = ps.evaluate_hessian(true_parameter, 1e-7, simulation_tuple)

        HMatrix[angleOneCounter, angleTwoCounter] = H
        detHMatrix[angleOneCounter, angleTwoCounter] = detH
        condHMatrix[angleOneCounter, angleTwoCounter] = condH
        detH0Matrix[angleOneCounter, angleTwoCounter] = detH0

        #destroy_routine(ps.simulation)
        angleTwoCounter += 1
        print angleTwoCounter

    angleOneCounter += 1
    angleTwoCounter = 0

# Once all the orientations have been analysed, plot the resulting Hessian metrics against the orientation design
# variables. This must be done in a separate visualisation script, so save the output data.
np.savetxt('Z:\opt\Visualisation\sensitivity_detH_output.txt', detHMatrix, delimiter=' || ', newline='\n')
designVariables = np.array([designVariableOneStart, designVariableOneFinish, designVariableOneStep, designVariableTwoStart, designVariableTwoFinish, designVariableTwoStep])
np.savetxt('Z:\opt\Visualisation\sensitivity_design_variables.txt', designVariables, delimiter=" || ", newline="\n")

print 'Done!'
