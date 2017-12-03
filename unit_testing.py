import unittest
import numpy as np
import CantileverSimulation

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
        self.assertTrue(function_which_gets_BC_type of node(1) == 'FIXED')
        self.assertTrue(function_which_gets_BC_type of node(3) == 'FIXED')
        self.assertTrue(function_which_gets_BC_type of node(5) == 'FIXED')
        self.assertTrue(function_which_gets_BC_type of node(7) == 'FIXED')



