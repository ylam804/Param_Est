from opencmiss.iron import iron
from parameter_optimisation import ParameterEstimation
from cantilever_simulation import CantileverSimulation
from cantilever_simulation import single_layer_objective_function
from cantilever_simulation import two_layer_objective_function
import numpy as np
import math

def destroy_routine(simulation):
    """
    Destroys some parts of the simulation so that they can be re-initialised with new values.
    """

    simulation.coordinate_system.Destroy()
    simulation.region.Destroy()
    simulation.basis.Destroy()
    simulation.problem.Destroy()

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


dimensions = np.array([30, 12, 12])
elements = np.array([2, 2, 2])
initial_parameter = np.array([8.4378])

theta1 = -90*math.pi/180
phi1 = 90*math.pi/180
theta2 = -90*math.pi/180
phi2 = 90*math.pi/180

ps = ParameterEstimation()
ps.simulation = CantileverSimulation()
ps.simulation.set_cantilever_dimensions(dimensions)
ps.simulation.set_cantilever_elements(elements)
ps.simulation.set_gravity_vector(ps.simulation.gravity_vector_calculation(theta1, phi1))
ps.simulation.set_diagnostic_level(1)
ps.simulation.setup_cantilever_simulation()
ps.simulation.set_Neo_Hookean_single_layer(initial_parameter)
ps.simulation.solve_simulation()

reset_problem(ps.simulation)
ps.simulation.set_Neo_Hookean_single_layer(initial_parameter+1)
ps.simulation.solve_simulation()

reset_problem(ps.simulation)
ps.simulation.set_Neo_Hookean_single_layer(initial_parameter+2)
ps.simulation.solve_simulation()
