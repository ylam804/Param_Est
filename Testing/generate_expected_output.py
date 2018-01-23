import numpy as np
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
from parameter_optimisation import ParameterEstimation

def generate_one_layer_cantilever_output():
    """
    Generates expected output data set for test_one_layer_cantilever_output unit test.

    Stiffness value is 5.432
    """

    # Set up the simulation
    sim = CantileverSimulation()
    sim.set_cantilever_dimensions(np.array([30, 12, 12]))
    sim.set_cantilever_elements(np.array([4, 4, 4]))
    sim.setup_cantilever_simulation()
    sim.set_Neo_Hookean_single_layer(np.array([5.432]))

    # Solve the simulation
    sim.solve_simulation()

    # Export the results
    sim.export_results("Expected Output/generate_one_layer_cantilever_output")

generate_one_layer_cantilever_output()
