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
import math
import numpy as np
from CantileverSimulation import CantileverSimulation

class ConvergenceTest:
    """
    Class for holding and running convergence tests of a FE model.
    """

    def __init__(self):
        """
        Create new instance of a test.
        """

        self.simulation =  None

    def set_simulation(self, simulation):
        """
        Sets the type of simulation model which will have the convergence tested on.

        :param simulation: A class containing all the set-up methods required to get the desired simulation running.
        :return:
        """

        self.simulation = simulation

    def destroy_routine(self):
        """
        Destroys some parts of the simulation so that they can be re-initialised with new values.
        """

        self.simulation.coordinate_system.Destroy()
        self.simulation.region.Destroy()
        self.simulation.basis.Destroy()
        self.simulation.problem.Destroy()

def destroy_routine(simulation):
    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()

def error_calculation(currentDataSet, previousDataSet):
    xError = yError = zError = np.zeros(len(currentDataSet))

    for i in range(len(currentDataSet)):
        xError[i] = abs(currentDataSet[i,0] - previousDataSet[i,0])
        yError[i] = abs(currentDataSet[i,1] - previousDataSet[i,1])
        zError[i] = abs(currentDataSet[i,2] - previousDataSet[i,2])

    xError = np.max(xError)
    yError = np.max(yError)
    zError = np.max(zError)

    RMSError = np.zeros(len(currentDataSet))
    for i in range(len(currentDataSet)):
        RMSError[i] = np.sqrt(np.average(np.array([(currentDataSet[i,0] - previousDataSet[i,0]) ** 2, (currentDataSet[i,1] - previousDataSet[i,1]) ** 2, (currentDataSet[i,2] - previousDataSet[i,2]) ** 2])))
    maxRMSError = np.max(RMSError)

    return xError, yError, zError, maxRMSError

def displacement_calculation(currentDataSet):
    displacementArray = np.zeros((1,len(currentDataSet)))

    for i in range(len(currentDataSet)):
        displacementArray[0,i] = math.sqrt(currentDataSet[i,0] ** 2 + currentDataSet[i,1] ** 2 + currentDataSet[i,2] ** 2)

    return displacementArray

# Set the tolerance required for the mesh convergence study
tolerance = 0.3

# Prepare the initial conditions for the first simulation
cantilever_dimensions = np.array([30, 12, 12])
cantilever_elements = np.array([1, 1, 1])
cantilever_initial_parameters = np.array([1.452])
simulation = CantileverSimulation()

# Now run the first simulation and collect the first data set which can then be used inside the while loop to
# determine if mesh convergence has been reached or not.
simulation.set_cantilever_dimensions(cantilever_dimensions)
simulation.set_cantilever_elements(cantilever_elements)
simulation.set_diagnostic_level(1)
simulation.setup_cantilever_simulation()
simulation.set_Mooney_Rivlin_parameter_values(cantilever_initial_parameters)

# Solve the simulation and derive the data from the results.
simulation.solve_simulation()
previousDataPoints = simulation.generate_data(3)
displacements = displacement_calculation(previousDataPoints)
print previousDataPoints
print '\n'

# A second data set is required to start the loop, so repeat another simulation after changing the number of elements
# only slightly. This requires clearing the previous simulation.
cantilever_elements = np.array([2, 2, 2])
destroy_routine(simulation)

simulation = CantileverSimulation()
simulation.set_cantilever_dimensions(cantilever_dimensions)
simulation.set_cantilever_elements(cantilever_elements)
simulation.set_diagnostic_level(1)
simulation.setup_cantilever_simulation()
simulation.set_Mooney_Rivlin_parameter_values(cantilever_initial_parameters)
simulation.solve_simulation()
currentDataPoints = simulation.generate_data(3)
displacements = np.append(displacements, displacement_calculation(currentDataPoints), axis = 0)
print currentDataPoints

# Calculate the error in the data in each of the axial directions using the two sets of data.
[xError, yError, zError, RMSError] = error_calculation(currentDataPoints, previousDataPoints)
xErrorArray = np.array([xError])
yErrorArray = np.array([yError])
zErrorArray = np.array([zError])
RMSErrorArray = np.array([RMSError])
# Now prepare the loop variables.
iteration = 0
converged = False

while converged == False and iteration < 10:
    # First, clear the previous simulation out of the simulation variable so it can be used again.
    destroy_routine(simulation)
    simulation = CantileverSimulation()

    # Now calculate the new variables which should be used to initialise the next simulation using the error in each
    # direction from the previous iteration.
    if xError > tolerance:
        cantilever_elements[0] = round(cantilever_elements[0] * 1.6)
    if yError > tolerance:
        cantilever_elements[1] = round(cantilever_elements[1] * 1.3)
    if zError > tolerance:
        cantilever_elements[2] = round(cantilever_elements[2] * 1.3)

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
    print '\n\n'
    iteration += 1
    print "Mesh Refinement = %d" % iteration
    print '\n'
    print previousDataPoints
    print '\n'
    print currentDataPoints

    # Calculate the overall error which will be used to check if the overall tolerance has been satisfied.
    errorArray = np.zeros((len(currentDataPoints), len(currentDataPoints[0])))
    errorCount = 0
    for i in range(len(currentDataPoints)):
        for j in range(len(currentDataPoints[i])):
            errorArray[i,j] = abs((currentDataPoints[i,j] - previousDataPoints[i,j]))
            if errorArray[i,j] > tolerance:
                errorCount += 1

    if errorCount == 0:
        converged = True

    displacements = np.append(displacements, displacement_calculation(currentDataPoints), axis = 0)

    # Calculate the error in the x, y and z directions.
    [xError, yError, zError, RMSError] = error_calculation(currentDataPoints, previousDataPoints)
    xErrorArray = np.append(xErrorArray, xError)
    yErrorArray = np.append(yErrorArray, yError)
    zErrorArray = np.append(zErrorArray, zError)
    RMSErrorArray = np.append(RMSErrorArray, RMSError)
    simulation.export_results()

# Return the final element dimensions required to converge to the tolerance
print '\n\n\n'
print "Final Number of Elements:"
print "              x = %d" % cantilever_elements[0]
print "              y = %d" % cantilever_elements[1]
print "              z = %d" % cantilever_elements[2]

print displacements

# Now construct an array which will be saved to a file for visualisation.
printingArray = np.array([xErrorArray])
printingArray = np.append(printingArray, np.array([yErrorArray]), axis=0)
printingArray = np.append(printingArray, np.array([zErrorArray]), axis=0)
printingArray = np.append(printingArray, np.array([RMSErrorArray]), axis=0)
np.savetxt('convergence_error_output.txt', printingArray, delimiter=' || ', newline="\n")
