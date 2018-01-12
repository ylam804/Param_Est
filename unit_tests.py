import unittest
import numpy as np
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import cantilever_objective_function

class TestCantileverSimulation(unittest.TestCase):

    def test_boundary_elements(self):
        """
        Tests that the Cantilever simulation is creating the boundary conditions of a simple cube as expected.
        """

        sim = CantileverSimulation()
        sim.set_cantilever_dimensions(np.array([30, 12, 12]))
        sim.set_cantilever_elements(np.array([8, 6, 8]))
        sim.setup_cantilever_simulation()

        # Now test the four nodes on the y-face of the element and check that they are all under the fixed condition.
        #self.assertTrue(function_which_gets_BC_type of node(1) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(3) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(5) == 'FIXED')
        #self.assertTrue(function_which_gets_BC_type of node(7) == 'FIXED')

