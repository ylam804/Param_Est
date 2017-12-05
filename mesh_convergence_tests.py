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



#import opencmiss.iron as iron
import cmath
import numpy as np
from CantileverSimulation import CantileverSimulation

def destroy_routine(simulation):
    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()

def error_calculation(currentDataSet, previousDataSet):
    xError = yError = zError = 0

    for i in range(len(currentDataSet)):
        xError += (currentDataSet[i,0] - previousDataSet[i,0]) ** 2
        yError += (currentDataSet[i,1] - previousDataSet[i,1]) ** 2
        zError += (currentDataSet[i,2] - previousDataSet[i,2]) ** 2

    xError = cmath.sqrt(xError/len(currentDataSet))
    yError = cmath.sqrt(yError/len(currentDataSet))
    zError = cmath.sqrt(zError/len(currentDataSet))

    return xError, yError, zError

# Set the tolerance required for the mesh convergence study
tolerance = 0.00001

# Prepare the initial conditions for the first simulation
cantilever_dimensions = np.array([60, 40, 40])
cantilever_elements = np.array([1, 1, 1])
cantilever_initial_parameters = np.array([1.5, 1.0])
simulation = CantileverSimulation()

# Now run the first simulation and collect the first data set which can then be used inside the while loop to
# determine if mesh convergence has been reached or not.
simulation.set_cantilever_dimensions(cantilever_dimensions)
simulation.set_cantilever_elements(cantilever_elements)
simulation.set_diagnostic_level(0)
simulation.setup_cantilever_simulation()
simulation.set_Mooney_Rivlin_parameter_values(cantilever_initial_parameters)

# Solve the simulation and derive the data from the results.
simulation.solve_simulation()
previousDataPoints = simulation.generate_data(3)
print previousDataPoints

# A second data set is required to start the loop, so repeat another simulation after changing the number of elements
# only slightly. This requires clearing the previous simulation.
cantilever_elements = np.array([2, 2, 2])
destroy_routine(simulation)

simulation = CantileverSimulation()
simulation.set_cantilever_dimensions(cantilever_dimensions)
simulation.set_cantilever_elements(cantilever_elements)
simulation.set_diagnostic_level(0)
simulation.setup_cantilever_simulation()
simulation.set_Mooney_Rivlin_parameter_values(cantilever_initial_parameters)
simulation.solve_simulation()
currentDataPoints = simulation.generate_data(3)
print currentDataPoints

# Calculate the error in the data in each of the axial directions using the two sets of data.
[xError, yError, zError] = error_calculation(currentDataPoints, previousDataPoints)

# Now prepare the loop variables.
iteration = 1
converged = False

while converged == False:
    # First, clear the previous simulation out of the simulation variable so it can be used again.
    destroy_routine(simulation)
    simulation = CantileverSimulation()

    # Now calculate the new variables which should be used to initialise the next simulation using the error in each
    # direction from the previous iteration.
    if xError > tolerance:
        cantilever_elements[0] = cantilever_elements[0] * round(cmath.log10(xError/tolerance))
    if yError > tolerance:
        cantilever_elements[1] = cantilever_elements[1] * round(cmath.log10(yError/tolerance))
    if yError > tolerance:
        cantilever_elements[2] = cantilever_elements[2] * round(cmath.log10(zError/tolerance))

    # Now set up and solve the next simulation with these parameters.
    simulation = CantileverSimulation()
    simulation.set_cantilever_dimensions(cantilever_dimensions)
    simulation.set_cantilever_elements(cantilever_elements)
    simulation.set_diagnostic_level(0)
    simulation.setup_cantilever_simulation()
    simulation.set_Mooney_Rivlin_parameter_values(cantilever_initial_parameters)
    simulation.solve_simulation()

    # Move the data from the previous simulation from the current to the previous variable before retrieving the
    # next set of data from the complete simulation.
    previousDataPoints = currentDataPoints
    currentDataPoints = simulation.generate_data(3)

    # Calculate the overall error which will be used to check if the overall tolerance has been satisfied.
    Error = 0
    for i in range(len(currentDataPoints)):
        for j in range(len(currentDataPoints[i])):
            Error += (currentDataPoints[i,j] - previousDataPoints[i,j]) ** 2

    Error = cmath.sqrt(zError/len(currentDataPoints))

    if Error < tolerance:
        converged = True

    # Calculate the error in the x, y and z directions.
    [xError, yError, zError] = error_calculation(currentDataPoints, previousDataPoints)

# Return the final element dimensions required to converge to the tolerance
print cantilever_elements
