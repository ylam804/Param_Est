import numpy as np
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
from scipy.optimize import least_squares


class ParameterEstimation:
    """
    Parameter estimation class
    """

    def __init__(self):
        """
        Create a new material law instance with no defaults.
        """
        self.lower_bounds = -np.inf
        self.upper_bounds = np.inf
        self.initial_parameters = None
        self.objective_function = None
        self.objective_function_arguments = ()
        self.solutions = None
        self.simulation = None

    def set_upper_bounds(self,upper_bounds=np.inf):
        """
        Set upper bounds. Defaults to upper_bounds=np.inf ie not bounds.
        Bounds need to be specified as a numpy 1D array with a length equal
        to the number of parameters being identified.
        """
        self.upper_bounds = upper_bounds

    def set_lower_bounds(self,lower_bounds=-np.inf):
        """
        Set lower bounds. Defaults to lower_bounds=np.inf ie not bounds.
        Bounds need to be specified as a numpy 1D array with a length equal
        to the number of parameters being identified.
        """
        self.lower_bounds = lower_bounds

    def set_initial_parameters(self,initial_parameters):
        """
        Set initial parameters (numpy array)
        """
        self.initial_parameters = initial_parameters

    def set_objective_function(self, objective_function, arguments=None):
        """
        Set the objective function
        """
        if callable(objective_function):
            self.objective_function = objective_function
        else:
            raise TypeError('The objective function must be callable ie it '
                            'must be a function')
        if arguments is not None:
            if isinstance(arguments, tuple):
                self.objective_function_arguments = arguments
            else:
                raise TypeError('The objective function arguments must be a '
                                'tuple')

    def optimise(self):
        """
        Routine for running the optimisation
        """
        if self.initial_parameters is None:
            raise ValueError('Initial parameters need to be defined')
        if self.objective_function is None:
            raise ValueError('An objective function need to be set')

        self.solutions = least_squares(self.objective_function,
                                       self.initial_parameters,
                                       args=self.objective_function_arguments,
                                       bounds=(self.lower_bounds,
                                               self.upper_bounds),
                                       diff_step=1e-5)
        return self.solutions

    def evaluate_hessian(self, x, stepsize):
        """
        Routine for evaluating the Hessian matrix using central finite differences
        """
        objfun = self.objective_function
        n = len(x)
        A = np.zeros(n)
        B = np.zeros(n)
        ee = stepsize * np.eye(n)

        # First-order derivatives: 2n function calls needed
        for i in range(n):
            A[i] = objfun(x + ee[:, i], self.simulation)
            B[i] = objfun(x - ee[:, i], self.simulation)

        # Second-order derivatives based on function calls only (Abramowitz and Stegun 1972, p.884): for dense Hessian, 2n+4n^2/2 function calls needed.
        H = np.zeros((n, n))
        for i in range(n):
            C = objfun(x + 2 * ee[:, i], self.simulation)
            E = objfun(x, self.simulation)
            F = objfun(x - 2 * ee[:, i], self.simulation)
            H[i, i] = (- C + 16 * A[i] - 30 * E + 16 * B[i] - F) / (12 * (ee[i, i] ** 2))
            for j in range(i + 1, n):
                G = objfun(x + ee[:, i] + ee[:, j], self.simulation)
                I = objfun(x + ee[:, i] - ee[:, j], self.simulation)
                J = objfun(x - ee[:, i] + ee[:, j], self.simulation)
                K = objfun(x - ee[:, i] - ee[:, j], self.simulation)
                H[i, j] = (G - I - J + K) / (4 * ee[i, i] * ee[j, j])
                H[j, i] = H[i, j]

        import cmath
        n = len(H)
        detH = np.linalg.det(H)
        condH = 1.0 / np.linalg.cond(H)
        H0 = np.zeros((n, n))
        for j in range(n):
            for k in range(n):
                H0[j, k] = H[j, k] / np.abs((cmath.sqrt(H[j, j] * H[k, k])))
        detH0 = np.linalg.det(H0)
        return H, detH, condH, detH0

###########
# Testing #
###########

if __name__ == "__main__":

    #data = np.array(([[54.127, 0.724, -11.211], [54.127, 39.276, -11.211], [64.432, -0.669, 27.737], [64.432, 40.669, 27.737]]))
    cantilever_dimensions = np.array([30, 12, 12])
    cantilever_elements = np.array([2, 2, 2])
    cantilever_initial_parameter = np.array([14.05])
    cantilever_guess_parameter = np.array([10.0])

    print "Parameters for data generation: C10 = 2.05, C01 = 0.0"
    print "Initial parameters for optimisation: C10 = 0.82, C01 = 0.0"

    ps = ParameterEstimation()
    ps.simulation = CantileverSimulation()
    ps.simulation.set_cantilever_dimensions(cantilever_dimensions)
    ps.simulation.set_cantilever_elements(cantilever_elements)
    ps.simulation.set_diagnostic_level(1)
    ps.simulation.setup_cantilever_simulation()
    ps.simulation.set_Neo_Hookean_single_layer(cantilever_initial_parameter)
    ps.simulation.solve_simulation()
    data = ps.simulation.generate_data(1)

    ps.simulation.set_projection_data(data)
    ps.initial_parameters = cantilever_guess_parameter
    simulation_tuple = (ps.simulation,)
    ps.set_objective_function(single_layer_objective_function, simulation_tuple)
    ps.optimise()
    print '\n'
    print "Optimisation routine results:"
    print ps.solutions.x

    [H, detH, condH, detH0] = ps.evaluate_hessian(ps.solutions.x, 1.e-7)
    print '\n'
    print "Hessian determinant"
    print detH

    #print  ps.objective_function(cantilever_initial_parameter, ps.simulation)

    #[H, detH, condH, detH0] = ps.evaluate_hessian(ps.solutions.x, 1e-7)
    #print detH
