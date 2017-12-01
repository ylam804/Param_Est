import numpy as np
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
            A[i] = objfun(x + ee[:, i])
            B[i] = objfun(x - ee[:, i])

        # Second-order derivatives based on function calls only (Abramowitz and Stegun 1972, p.884): for dense Hessian, 2n+4n^2/2 function calls needed.
        H = np.zeros((n, n))
        for i in range(n):
            C = objfun(x + 2 * ee[:, i])
            E = objfun(x)
            F = objfun(x - 2 * ee[:, i])
            H[i, i] = (- C + 16 * A[i] - 30 * E + 16 * B[i] - F) / (12 * (ee[i, i] ** 2))
            for j in range(i + 1, n):
                G = objfun(x + ee[:, i] + ee[:, j])
                I = objfun(x + ee[:, i] - ee[:, j])
                J = objfun(x - ee[:, i] + ee[:, j])
                K = objfun(x - ee[:, i] - ee[:, j])
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


def testParameterEstimation():

    """
    Test the parameter estimation routines
    """

    ps = ParameterEstimation()
    ps.set_initial_parameters(np.array([0.5,0.5]))
    ps.set_objective_function(fun_rosenbrock_mse())
    ps.optimise()
    ps.set_objective_function(fun_rosenbrock_mse)
    ps.H, ps.detH, ps.condH, ps.detH0 = ps.evaluate_hessian(ps.solutions.x, 1.e-7)

    return ps


def fun_rosenbrock(x):
    """
    An objective function for testing
    """
    return np.array([10 * (x[1] - x[0] ** 2), (1 - x[0])])


def fun_rosenbrock_mse(x):
    """
    An objective function for testing
    """
    vector = np.array([100 * (x[1] - x[0] ** 2), (1 - x[0])])
    mse = np.mean(vector * vector)
    return mse


def plotParameterEstimationError(x, y, dx, dy):
    """
    Using a grid of points as some initial inputs to the parameter estimation function, plots the error at each of
    points.

    :param x: a np.array containing two numbers which are the initial and final values of the first parameter.
    :param y: a np.array containing two numbers which are the initial and final values of the second parameter.
    :param dx: the number of points needed in the x direction.
    :param dy: the number of points needed in the y direction.
    :return: None, but plots a graph of the error as the two parameters change
    """

    # First create the grid space by making an array which extends from the initial to final value of each parameter
    # in steps of the grid spacing.
    xvalues = np.linspace(x[0], x[1], dx)
    yvalues = np.linspace(y[0], y[1], dy)

    # Create a matrix which can store all the error values calculated by the ParameterEstimation function.
    errorMatrix = np.zeros((dx, dy))

    # Now loop through these two arrays and call the ParameterEstimation function to find the error associated with
    # each combination of initial parameter values.
    for i in range(dx):
        for j in range(dy):
            ps = ParameterEstimation()
            ps.set_initial_parameters(np.array([xvalues[i],yvalues[j]]))
            ps.set_objective_function(fun_rosenbrock_mse)
            errorMatrix[i,j] = ps.objective_function(ps.initial_parameters)



    # Now that the error has been calculated at each point in the grid, plot this as a surface.
    #import matplotlib.pyplot as plt
    #from matplotlib import cm
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #xgrid, ygrid = np.meshgrid(yvalues, xvalues)
    #surf = ax.plot_surface(xgrid, ygrid, errorMatrix, cmap=cm.coolwarm,
    #                       linewidth=0, antialiased=False)
    #plt.show()

    # Now calculate the optimal values of
    ps.optimise()
    ps.set_objective_function(fun_rosenbrock_mse)
    ps.H, ps.detH, ps.condH, ps.detH0 = ps.evaluate_hessian(ps.solutions.x, 1.e-7)

###########
# Testing #
###########


ps = ParameterEstimation
s.set_initial_parameters(np.array([0.5,0.5]))
ps.set_objective_function(cantilever_objective_function)
ps.optimise()

