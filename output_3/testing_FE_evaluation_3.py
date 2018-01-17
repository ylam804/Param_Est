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

    #simulation.controlLoop.Destroy()
    #simulation.boundaryConditions.Destroy()
    #simulation.decomposition.Destroy()
    #simulation.equations.Destroy()
    #simulation.generatedMesh.Destroy()
    #simulation.geometricField.Destroy()
    #simulation.fibreField.Destroy()
    #simulation.dependentField.Destroy()
    #simulation.materialField.Destroy()
    #simulation.sourceField.Destroy()
    #simulation.mesh.Destroy()
    #simulation.coordinate_system.Destroy()
    #simulation.basis.Destroy()
    #simulation.problem.Destroy()
    #simulation.region.Destroy()


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

#ps.simulation.export_results("output_3/Cantilever")
ps.simulation.set_projection_data()
ps.set_objective_function(single_layer_objective_function)
[H, detH, condH, detH0] = ps.evaluate_hessian(initial_parameter, 1e-7)

print "For angles Theta = {0}, Phi = {1}".format(theta1, phi1)
print "Gravity Vector = {0}".format(ps.simulation.gravity_vector)
print "Determinant = {0}".format(detH)
print "\n"


destroy_routine(ps.simulation)
ps = None

ps = ParameterEstimation()
ps.simulation = CantileverSimulation()
ps.simulation.set_cantilever_dimensions(dimensions)
ps.simulation.set_cantilever_elements(elements)
ps.simulation.set_gravity_vector(ps.simulation.gravity_vector_calculation(theta2, phi2))
ps.simulation.set_diagnostic_level(1)
ps.simulation.setup_cantilever_simulation()
ps.simulation.set_Neo_Hookean_single_layer(initial_parameter)
ps.simulation.solve_simulation()

#ps.simulation.export_results("output_3/Cantilever")
ps.simulation.set_projection_data()
ps.set_objective_function(single_layer_objective_function)
[H, detH, condH, detH0] = ps.evaluate_hessian(initial_parameter, 1e-7)

print "For angles Theta = {0}, Phi = {1}".format(theta2, phi2)
print "Gravity Vector = {0}".format(ps.simulation.gravity_vector)
print "Determinant = {0}".format(detH)

destroy_routine(ps.simulation)
ps = None
