from opencmiss.iron import iron
from cantilever_simulation import CantileverSimulation
from parameter_optimisation import ParameterEstimation
import numpy as np
import cmath


def reset_problem(simulation):
    """
    Resets the FE model's solvers so there is no data from the previous solve used in the next one.

    :param simulation: The simulation which needs its solver reset.
    """

    problemUserNumber = 12
    numberOfLoadIncrements = 1
    NumberOfGaussXi = 4

    simulation.problem.Destroy()

    simulation.problem = iron.Problem()
    simulation.problemSpecification = [iron.ProblemClasses.ELASTICITY,
            iron.ProblemTypes.FINITE_ELASTICITY,
            iron.ProblemSubtypes.NONE]
    simulation.problem.CreateStart(problemUserNumber, simulation.problemSpecification)
    simulation.problem.CreateFinish()

    # Create the problem control loop
    simulation.problem.ControlLoopCreateStart()
    simulation.controlLoop = iron.ControlLoop()
    simulation.problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],simulation.controlLoop)
    simulation.controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
    simulation.problem.ControlLoopCreateFinish()

    simulation.nonLinearSolver = iron.Solver()
    simulation.linearSolver = iron.Solver()
    simulation.problem.SolversCreateStart()
    simulation.problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,simulation.nonLinearSolver)

    if simulation.diagnostics == 4 or simulation.diagnostics == 'Matrix':
        simulation.nonLinearSolver.outputType = iron.SolverOutputTypes.MATRIX
    elif simulation.diagnostics == 3 or simulation.diagnostics == 'Solver':
        simulation.nonLinearSolver.outputType = iron.SolverOutputTypes.SOLVER
    elif simulation.diagnostics == 2 or simulation.diagnostics == 'Timing':
        simulation.nonLinearSolver.outputType = iron.SolverOutputTypes.TIMING
    elif simulation.diagnostics == 1 or simulation.diagnostics == 'Progress':
        simulation.nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
    else:
        simulation.nonLinearSolver.outputType = iron.SolverOutputTypes.NONE

    simulation.nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
    simulation.nonLinearSolver.NewtonAbsoluteToleranceSet(1e-9)
    simulation.nonLinearSolver.NewtonSolutionToleranceSet(1e-9)
    simulation.nonLinearSolver.NewtonRelativeToleranceSet(1e-9)
    simulation.nonLinearSolver.NewtonMaximumIterationsSet(int(1e6))
    simulation.nonLinearSolver.NewtonMaximumFunctionEvaluationsSet(int(1e6))
    simulation.nonLinearSolver.NewtonLinearSolverGet(simulation.linearSolver)
    simulation.linearSolver.linearType = iron.LinearSolverTypes.DIRECT
    #linearSolver.libraryType = iron.SolverLibraries.LAPACK
    simulation.problem.SolversCreateFinish()

    simulation.solver = iron.Solver()
    simulation.solverEquations = iron.SolverEquations()
    simulation.problem.SolverEquationsCreateStart()
    simulation.problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,simulation.solver)
    simulation.solver.SolverEquationsGet(simulation.solverEquations)
    simulation.solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    equationsSetIndex = simulation.solverEquations.EquationsSetAdd(simulation.equationsSet)
    simulation.problem.SolverEquationsCreateFinish()

    simulation.boundaryConditions = iron.BoundaryConditions()
    simulation.solverEquations.BoundaryConditionsCreateStart(simulation.boundaryConditions)

    numberOfNodes = (NumberOfGaussXi + (NumberOfGaussXi-1)*(simulation.cantilever_elements[0]-1))\
                    * (NumberOfGaussXi + (NumberOfGaussXi-1)*(simulation.cantilever_elements[1]-1))\
                     * (NumberOfGaussXi + (NumberOfGaussXi-1)*(simulation.cantilever_elements[2]-1))

    for nodeNum in range(1, numberOfNodes+1, (NumberOfGaussXi + (NumberOfGaussXi-1)*(simulation.cantilever_elements[0]-1))):
        simulation.boundaryConditions.AddNode(simulation.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        simulation.boundaryConditions.AddNode(simulation.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
        simulation.boundaryConditions.AddNode(simulation.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

    simulation.solverEquations.BoundaryConditionsCreateFinish()

def reset_deformed_field(simulation):
    """

    :param simulation:
    :return:
    """

    dependentFieldUserNumber = 13

    simulation.dependentField.Destroy()

    simulation.dependentField = iron.Field()
    simulation.equationsSet.DependentCreateStart(dependentFieldUserNumber,simulation.dependentField)
    simulation.dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
    simulation.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
    simulation.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
    #if UsePressureBasis:
    #    # Set the pressure to be nodally based and use the second mesh component
    #    if InterpolationType == 4:
    #        simulation.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
    #        simulation.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
    #    simulation.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
    #    simulation.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
    #if InterpolationType == 4:
    #    simulation.dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    simulation.equationsSet.DependentCreateFinish()

    iron.Field.ParametersToFieldParametersComponentCopy(
        simulation.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        simulation.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        simulation.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        simulation.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        simulation.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        simulation.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    iron.Field.ComponentValuesInitialiseDP(
        simulation.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0)

def create_routine(simulation, dimensions, elements, gravity_vector, diagnostic_level, data):
    """

    :param dimensions:
    :param elements:
    :param gravity_vector:
    :param diagnostic_level:
    :param data:
    :return:
    """

    destroy_routine(simulation)
    simulation = None

    simulation = CantileverSimulation()
    simulation.set_cantilever_dimensions(dimensions)
    simulation.set_cantilever_elements(elements)
    simulation.set_gravity_vector(gravity_vector)
    simulation.set_diagnostic_level(diagnostic_level)
    simulation.setup_cantilever_simulation()
    simulation.data = data

    return simulation

def destroy_routine(simulation):
    """
    Destroys some parts of the simulation so that they can be re-initialised with new values.
    """

    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()

def new_evaluate_hessian_method(objective_function, x, stepSize, simulation):
    """
    Routine for evaluating the Hessian matrix using central finite differences
    """

    objfun = objective_function
    n = len(x)
    A = np.zeros(n)
    B = np.zeros(n)
    ee = stepSize * np.eye(n)

    # First-order derivatives: 2n function calls needed
    for i in range(n):
        reset_problem(simulation)
        reset_deformed_field(simulation)
        A[i] = objfun(x + ee[:, i], simulation)
        B[i] = objfun(x - ee[:, i], simulation)

    # Second-order derivatives based on function calls only (Abramowitz and Stegun 1972, p.884): for dense Hessian, 2n+4n^2/2 function calls needed.
    H = np.zeros((n, n))
    for i in range(n):
        C = objfun(x + 2 * ee[:, i], simulation)
        E = objfun(x, simulation)
        F = objfun(x - 2 * ee[:, i], simulation)
        H[i, i] = (- C + 16 * A[i] - 30 * E + 16 * B[i] - F) / (12 * (ee[i, i] ** 2))
        for j in range(i + 1, n):
            G = objfun(x + ee[:, i] + ee[:, j], simulation)
            I = objfun(x + ee[:, i] - ee[:, j], simulation)
            J = objfun(x - ee[:, i] + ee[:, j], simulation)
            K = objfun(x - ee[:, i] - ee[:, j], simulation)
            H[i, j] = (G - I - J + K) / (4 * ee[i, i] * ee[j, j])
            H[j, i] = H[i, j]

    n = len(H)
    detH = np.linalg.det(H)
    condH = 1.0 / np.linalg.cond(H)
    H0 = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            H0[j, k] = H[j, k] / np.abs((cmath.sqrt(H[j, j] * H[k, k])))
    detH0 = np.linalg.det(H0)
    return H, detH, condH, detH0

