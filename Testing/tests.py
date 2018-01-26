import unittest
from opencmiss.iron import iron
import numpy as np
import math
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
from parameter_optimisation import ParameterEstimation

def destroy_routine(simulation):
    """
    Destroy's a simulation which has been set up so that the next test can create its own simulation.

    :param simulation: The simulation to be destroyed.
    """

    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.pressureBasis.Destroy()
    simulation.problem.Destroy()

def test_set_diagnostic_level():
    """ Test if the set_diagnostic_level method sets the simulation's diagnostic level as expected. """

    print "\n"
    print "Testing Cantilever Simulation set_diagnostics_level function:"
    print "\n"

    # Create the simulation
    sim = CantileverSimulation()

    # Check its initial diagnostic level
    if (sim.diagnostics == None):
        print "Diagnostics Initiation Successful."
    else:
        print "Diagnostics Initiation Failed."

    # Call the set_diagnostic_level function and check if it changed the level correctly.
    sim.set_diagnostic_level(1)
    if (sim.diagnostics == 1):
        print "set_diagnostic_level function Successful."
    else:
        print "set_diagnostic_level function Failed."

    print "\n"


def test_gravity_vector_calculation():
    """ Tests that the function correctly calculates the direction of the gravity vector from two angles. """

    print "\n"
    print "Testing Cantilever Simulation gravity_vector_calculation function:"
    print "\n"

    # Create the simulation
    sim = CantileverSimulation()
    passCounter = 0

    # Run the tests
    if (sim.gravity_vector_calculation(90*math.pi/180, 0*math.pi/180).all() ==  np.array([-9.81, 0.0, 0.0]).all()):
        passCounter += 1
    if (sim.gravity_vector_calculation(90*math.pi/180, 0*math.pi/180).all() ==  np.array([-9.81, 0.0, 0.0]).all()):
        passCounter += 1

    # Can easily add more test cases!

    print "Calculation passed {0} of 2 tests.".format(passCounter)
    print "\n"

def test_node_moving_calculation():
    """ Tests that the function which moves the nodes of the FE model does so in the correct manner. """

    print "\n"
    print "Testing Cantilever Simulation node relocation function:"
    print "\n"

    # Create the simulation
    sim = CantileverSimulation()

    # Set up the counter
    passCounter = 0

    # Run the tests
    sim.set_cantilever_dimensions(np.array([30, 12, 12]))
    if (np.isclose(sim.calculate_new_x_value(30.0), 30)):
        passCounter += 1
    if (np.isclose(sim.calculate_new_x_value(0.0), 0)):
        passCounter += 1
    if (np.isclose(sim.calculate_new_x_value(15.0), 7.5)):
        passCounter += 1
    if (np.isclose(sim.calculate_new_x_value(20.0), 120.0/9.0)):
        passCounter += 1

    print "Calculation passed {0} of 4 tests.".format(passCounter)
    print "\n"

def test_node_moving_function():
    """ Tests whether the node moving function successfully moves the nodes closer to the model's fixed end. """

    print "\n"
    print "Testing Cantilever Simulation node relocation function:"
    print "\n"

    # Set up the simulation, making it small so the location of nodes can be easily determined.
    sim = CantileverSimulation()
    sim.set_cantilever_dimensions(np.array([30, 12, 12]))
    sim.set_cantilever_elements(np.array([1, 1, 1]))
    sim.setup_cantilever_simulation()

    # Since the simulation uses Cubic Lagrange interpolation, there are 64 nodes in the one element.
    sim.move_one_layer_nodes()

    # We know that the nodes should've been at 0, 10, 20 and 30 mm along the x-axis.
    passCounter = 0

    if np.isclose(sim.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, 1, 1), 0):
        passCounter += 1
    if np.isclose(sim.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, 2, 1), 30.0/9.0):
        passCounter += 1
    if np.isclose(sim.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, 3, 1), 120.0/9.0):
        passCounter += 1
    if np.isclose(sim.dependentField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,1, 1, 4, 1), 30):
        passCounter += 1

    destroy_routine(sim)

    print "Calculation passed {0} of 4 tests.".format(passCounter)
    print "\n"



test_set_diagnostic_level()
test_gravity_vector_calculation()
test_node_moving_calculation()
test_node_moving_function()
