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
import numpy as np
from CantileverSimulation import CantileverSimulation

def destroy_routine(simulation):
    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()



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

# Now prepare the loop variables.
iteration = 1
converged = False

while converged == False:
    # First, clear the previous simulation out of the simulation variable so it can be used again.
    simulation = None
    simulation = CantileverSimulation()

    # Now calculate the new variables which should be used to initialise the next simulation using the error in each
    # direction from the previous iteration.
