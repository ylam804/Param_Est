import unittest
import numpy as np
import CantileverSimulation
from CantileverSimulation import cantilever_objective_function

class CantileverBoundaryConditions(unittest.TestCase):

    def test_boundary_elements(self):
        """
        Tests that the Cantilever simulation is creating the boundary conditions of a simple cube as expected.
        """

        sim = CantileverSimulation()
        sim.set_cantilever_dimensions(np.array([60, 40, 40]))
        sim.set_cantilever_elements(np.array([1, 1, 1]))
        sim.setup_cantilever_simulation()

        # Now test the four nodes on the y-face of the element and check that they are all under the fixed condition.
        #self.assertTrue(function_which_gets_BC_type of node(1) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(3) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(5) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(7) == 'FIXED')


class CantileverDataGeneration(unittest.TestCase):

    def test_data_generation(self):
        """
        Tests that the four corners of the generated data are located at the corners when there is only one element in
        the FE model.
        """

        sim = CantileverSimulation()
        sim.set_cantilever_dimensions(np.array([60, 40, 40]))
        sim.set_cantilever_elements(np.array([1, 1, 1]))
        sim.set_gravity_vector(np.array([0.0, 0.0, 0.0]))
        sim.setup_cantilever_simulation()
        sim.prepare_projection()
        cantilever_objective_function(np.array([2.0, 0.0]), sim)
        dataLocations = sim.generate_data(3)

        self.assertTrue(dataLocations[0] == [60,0,0])
        self.assertTrue(dataLocations[1] == [60,40,0])
        self.assertTrue(dataLocations[2] == [60,0,40])
        self.assertTrue(dataLocations[3] == [60,40,40])


###########
# TESTING #
###########

dataTest = CantileverDataGeneration()
dataTest.test_data_generation()
