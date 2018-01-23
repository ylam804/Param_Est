import unittest
import numpy as np
import math
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
from parameter_optimisation import ParameterEstimation

def test_set_diagnostic_level():
    """
    Test if the set_diagnostic_level method sets the simulation's diagnostic level as expected.
    """

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
    """
    Tests that the function correctly calculates the direction of the gravity vector from two angles.
    """

    print "\n"
    print "Testing Cantilever Simulation gravity_vector_calculation function:"
    print "\n"

    # Create the simulation
    sim = CantileverSimulation()
    passCounter = 0

    # Run the tests
    if (sim.gravity_vector_calculation(90*math.pi/180, 0*math.pi/180)) == np.array([-9.81, 0.0, 0.0]):
        passCounter += 1
    if (sim.gravity_vector_calculation(0*math.pi/180, 90*math.pi/180)) == np.array([0.0, 9.81, 0.0]):
        passCounter += 1
    if (sim.gravity_vector_calculation(90*math.pi/180, 0*math.pi/180)) == np.array([-9.81, 0.0, 0.0]):
        passCounter += 1
    if (sim.gravity_vector_calculation(90*math.pi/180, 0*math.pi/180)) == np.array([-9.81, 0.0, 0.0]):
        passCounter += 1


test_set_diagnostic_level()
test_gravity_vector_calculation()
