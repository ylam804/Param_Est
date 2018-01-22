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

import math
import numpy as np
from cantilever_simulation import CantileverSimulation

class ConvergenceTest:
    """
    Class for holding and running convergence tests of a FE model.
    """

    def __init__(self):
        """
        Create new instance of a test.
        """

        self.sim =  None
        self.tolerance = 1e-3
        self.convergenceCriteria = False
        self.meshIterationCounter = 1
        self.elements = np.array([1, 1, 1])
        self.currentDataSet = None
        self.previousDataSet = None
        self.RMSError = self.tolerance * 2
        self.axialError = np.array([1.0, 1.0, 1.0])
        self.dataRecord = np.array([])
        self.elementRecord = np.array([])

    def set_simulation(self, simulation):
        """
        Sets the type of simulation model which will have the convergence tested on.

        :param simulation: A class containing all the set-up methods required to get the desired simulation running.
        """

        self.simulation = simulation

    def destroy_routine(self):
        """
        Destroys some parts of the simulation so that they can be re-initialised with new values.
        """

        self.sim.coordinate_system.Destroy()
        self.sim.region.Destroy()
        self.sim.basis.Destroy()
        self.sim.pressureBasis.Destroy()
        self.sim.problem.Destroy()

    def calculate_axial_error(self):
        """
        Calculates the error in each dimension at each point in the data set
        """

        axErr = np.zeros((3, len(self.currentDataSet)))

        for i in range(len(self.currentDataSet)):
            axErr[0, i] = ((self.currentDataSet[i,0] - self.previousDataSet[i,0])**2)
            axErr[1, i] = ((self.currentDataSet[i,1] - self.previousDataSet[i,1])**2)
            axErr[2, i] = ((self.currentDataSet[i,2] - self.previousDataSet[i,2])**2)

        self.axialError[0] = np.sqrt(np.average(axErr[0]))
        self.axialError[1] = np.sqrt(np.average(axErr[1]))
        self.axialError[2] = np.sqrt(np.average(axErr[2]))

    def calculate_RMS_error(self):
        """
        Find the RMS error of the absolute error distance between the current and previous corner pairs.
        """

        err = np.array([0.0, 0.0, 0.0, 0.0])

        for i in range(len(self.currentDataSet)):
            err[i] = np.sqrt((self.currentDataSet[i,0] - self.previousDataSet[i,0])**2
                             + (self.currentDataSet[i,1] - self.previousDataSet[i,1])**2
                             + (self.currentDataSet[i,2] - self.previousDataSet[i,2])**2)

        self.RMSError = np.sqrt(np.average(err))

    def calculate_new_elements(self, direction):
        """
        Calculate the new number of elements for a given direction.

        :param direction: The x, y, or z direction denoted by an integer of 0, 1, or 2 respectively.
        """

        newEls = self.elements[direction] * (self.axialError[direction]/self.tolerance) # Needs refining

        if newEls < self.elements[direction]:
            newEls = self.elements[direction]
        elif newEls > (self.elements[direction] + 2):
            newEls = self.elements[direction] + 1

        self.elements[direction] = newEls

    def store_data(self):
        """
        In order to calculate the errors, each iteration needs to be compared to the most recent one. To do this, every
        each iteration's data needs to be recorded. Then, once the test has ended each iteration can be compared and
        the residuals plotted.
        """

        if len(self.dataRecord) == 0:
            self.dataRecord = np.array([self.currentDataSet])
        else:
            self.dataRecord = np.append(self.dataRecord, np.array([self.currentDataSet]), axis=1)

    def store_elements(self):
        """
        Calculate and record the number of elements to be displayed on the plot of mesh convergence.
        """

        if len(self.elementRecord) == 0:
            self.elementRecord = np.array([self.elements[0] * self.elements[1] * self.elements[2]])
        else:
            self.elementRecord = np.append(self.elementRecord, [self.elements[0] * self.elements[1] * self.elements[2]])

# Create instance of ConvergenceTest class
conTest = ConvergenceTest()

# Define some useful variables.
dimensions = np.array([30, 12, 12])
parameterValue = np.array([7.589])
conTest.tolerance = 1e-3

# Add a simulation to the convergence
conTest.sim = CantileverSimulation()

# Set up the chosen simulation.
conTest.sim.set_cantilever_elements(np.array([2, 2, 2]))
conTest.sim.set_cantilever_dimensions(dimensions)
conTest.sim.set_diagnostic_level(1)
conTest.sim.setup_cantilever_simulation()
conTest.sim.set_Neo_Hookean_single_layer(parameterValue)

# Now solve the simulation and generate a data set
conTest.sim.solve_simulation()
conTest.currentDataSet = conTest.sim.generate_data(0)
conTest.store_data()
conTest.store_elements()

# Increase the number of elements before running the next simulation.
conTest.elements = np.array([3, 3, 3])

# Lastly create an array for storing the errors from each iteration so they can be plotted later.
errorArray = np.array([[1, 1, 1, 1]])

# Now start the convergence loop.
while conTest.meshIterationCounter < 6 and conTest.RMSError > conTest.tolerance:

    # First, reset the simulation.
    conTest.destroy_routine()
    conTest.sim = CantileverSimulation()

    # Now repeat all the settings.
    conTest.sim.set_cantilever_dimensions(np.array([30, 12, 12]))
    conTest.sim.set_cantilever_elements(conTest.elements)
    conTest.sim.set_diagnostic_level(1)
    conTest.sim.setup_cantilever_simulation()
    conTest.sim.set_Neo_Hookean_single_layer(np.array([7.589]))

    # Solve the simulation.
    conTest.sim.solve_simulation()

    # Move the previous iteration's generated data to the correct variable, then generate new data.
    conTest.previousDataSet = conTest.currentDataSet
    conTest.currentDataSet = conTest.sim.generate_data(0)
    conTest.store_data()
    conTest.store_elements()

    # Use these two data sets to calculate both the RMS error in each axial direction along with the RMS error for the
    # total distance between a corner in the two data sets.
    conTest.calculate_axial_error()
    conTest.calculate_RMS_error()

    # Using these errors, calculate how much the number of elements should increase by.
    for i in range(len(conTest.sim.cantilever_elements)):
        conTest.calculate_new_elements(i)

    # Then compile these errors into a 1-by-4 array and append them to the errorArray for printing to a file and
    # visualisation later.
    errorArray = np.append(errorArray, np.array([[conTest.RMSError, conTest.axialError[0], conTest.axialError[1], conTest.axialError[2]]]), axis=0)

    print "\n\n\n\n"
    print "Mesh Iteration {0} Complete".format(conTest.meshIterationCounter)
    print "\n\n\n\n"

    np.savetxt('convergence_data_record.txt', conTest.dataRecord[0], newline="\n")
    np.savetxt('convergence_element_record.txt', conTest.elementRecord, newline="\n")

    # Increase the convergence iteration counter.
    conTest.meshIterationCounter += 1

# Now that the convergence has finished, print out the array full of errors.
np.savetxt('convergence_error_output.txt', errorArray, newline="\n")

