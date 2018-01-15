from opencmiss.iron import iron
import numpy as np
import math
from input_output import exportDatapointsErrorExdata


class CantileverSimulation:
    """
    Class for simulating a cantilever.
    """

    def __init__(self):
        """
        Create a new simulation for the deformation of a cantilever beam.
        """

        self.data = None
        self.diagnostics = None
        self.initial_parameters = None
        self.cantilever_dimensions = None  # in mm
        self.cantilever_elements = None
        self.density = 9.0E-4  # in g mm^-3
        self.gravity_vector = [0.0, 0.0, -9.81]  # in m s^-2
        self.nodes = None
        self.coordinate_system = None
        self.region = None
        self.basis = None
        self.pressureBasis = None
        self.generatedMesh = None
        self.mesh = None
        self.decomposition = None
        self.geometricField = None
        self.fibreField = None
        self.equationsSetField = None
        self.equationsSet = None
        self.equationsSetSpecification = None
        self.dependentField = None
        self.materialField = None
        self.sourceField = None
        self.equations = None
        self.problem = None
        self.problemSpecification = None
        self.controlLoop = None
        self.nonLinearSolver = None
        self.linearSolver = None
        self.solver = None
        self.solverEquations = None
        self.boundaryConditions = None
        self.dataPoints = None
        self.dataProjection = None
        self.error = 0
        self.fields = None
        self.numPointsPerFace = 3

    def set_projection_data(self, data=None):
        """
        Initialises the set of data points which were scanned on the surface of the deformed experiment.

        :param data: A 2D array where each row is three point in mm of the x, y, and z coordinates respectively.
        """

        if data == None:
            self.data = self.generate_data(1)
        else:
            self.data = data

    def set_diagnostic_level(self, diagnostics_choice):
        """
        Set the level of diagnostic display required.

        :param diagnostics_choice: None = 0
                                    Progress = 1
                                    Timing = 2
                                    Solver = 3
                                    Matrix = 4
        """

        self.diagnostics = diagnostics_choice


    def set_initial_parameters(self, initial_parameters):
        """
        Set the initial cantilever parameters.

        :param initial_parameters: An array of values which will determine the values of the material parameters for
                the cantilever.
        """

        self.initial_parameters = initial_parameters

    def set_cantilever_dimensions(self, dimensions):
        """
        Set the overall dimensions of the cantilever beam.

        :param dimensions: An np array of three values, the first is the x-dimension, the second the y-dimension and
                            the third is the z-dimension.
        """

        self.cantilever_dimensions = dimensions

    def set_cantilever_elements(self, elements):
        """
        Set the number of elements in each dimension of the cantilever.

        :param elements: An np array of three values, the first is the number of elements in the x-dimension, the second
                            is the number in the the y-dimension and the third is the number in the z-dimension.
        """

        self.cantilever_elements = elements

    def set_gravity_vector(self, gravity_vector):
        """
        Define the direction of gravity for the FE simulation.

        :param gravity_vector: A 2D/3D vector with the magnitude gravity acting in each coordinate direction.
        """

        self.gravity_vector = gravity_vector

    def gravity_vector_calculation(self, theta, phi):
        """
        Given spherical coordinates, calculates the apparent direction of gravity experienced by the cantilever.

        :param theta: Elevation of the beam from the horizontal.
        :param phi: Twist of the beam about its longitudinal axis.
        """

        # To find the 3D rotation of the gravity vector, start by calculating the x and z components of the vector which
        # are due to the lift away from the x-axis denoted by theta.
        value = (-9.81) ** 2 * math.sin(theta) ** 2
        if theta > 0:
            v = np.array([-math.sqrt(value), 0, -(math.sqrt((-9.81) ** 2 - value))])
        else:
            v = np.array([math.sqrt(value), 0, -(math.sqrt((-9.81) ** 2 - value))])

        # Now, to find how this is influenced by the rotation of the beam about its longitudinal axis, use Rodrigues'
        # rotation formula to decompose that vector into its parts parallel and perpendicular to the beam and then only
        # rotate those parts which are perpendicular.

        # The formula is:
        # Vrotated = V cos(theta) + (k x V) * sin(theta) + k (k . V)(1 - cos(theta)

        # Where V is the unrotated vector.
        #       Vrotated is the resultant vector.
        #       k is a unit vector along the axis being rotated about.
        #       theta is the angle of rotation.

        # Since we know that the x-axis will be the one rotated about, we can set k
        k = np.array([1, 0, 0])

        # Calculate the cross of k and v
        kcrossV = np.array([(k[1]*v[2] - k[2]*v[1]), (k[2]*v[0] - k[0]*v[2]), (k[0]*v[1] - k[1]*v[0])])

        # Now calculate the dot of k and v
        kdotV = k[0]*v[0] + k[1]*v[1] + k[2]*v[2]

        # Now in the final calculation we replace the angel theta in the formula with phi because that is the angle
        # which denotes the rotation of the beam about its longitudinal axis.
        v = v * math.cos(phi) + (kcrossV * math.sin(phi)) + (k * kdotV * (1 - math.cos(phi)))

        return v

    def set_Xi_points_num(self, numPoints):
        """
        Sets the number of points along each side of a face which will be used to generate Xi coordinates used for
        interpolating the surface position of each element.

        This won't work with two or fewer points per side because one of the data projection created later would have
        zero points in it, which isn't allowed. So if the number of points is less than two, change it so that it is
        equal to three.

        :param numPointsPerFace: An integer value for the number of points along one side.
        """

        if numPoints <= 2:
            numPoints = 3

        self.numPointsPerFace = numPoints

    def set_cantilever_density(self, density):
        """
        Set the density of the cantilever beam.

        :param density: A float for the density in grams per millimeter cubed (g mm^-3)
        """

    def setup_cantilever_simulation(self):
        """
        Uses the input values above to prepare the cantilever simulation to be run.
        """

        # Set problem parameters
        width = self.cantilever_dimensions[0]
        length = self.cantilever_dimensions[1]
        height = self.cantilever_dimensions[2]
        density = self.density
        gravity = self.gravity_vector

        UsePressureBasis = False
        NumberOfGaussXi = 4

        coordinateSystemUserNumber = 1
        regionUserNumber = 1
        basisUserNumber = 1
        pressureBasisUserNumber = 2
        generatedMeshUserNumber = 1
        meshUserNumber = 1
        decompositionUserNumber = 1
        geometricFieldUserNumber = 1
        fibreFieldUserNumber = 2
        materialFieldUserNumber = 3
        dependentFieldUserNumber = 4
        sourceFieldUserNumber = 5
        equationsSetFieldUserNumber = 6
        equationsSetUserNumber = 1
        problemUserNumber = 1

        # Set all diagnostic levels on for testing
        # iron.DiagnosticsSetOn(iron.DiagnosticTypes.All,[1,2,3,4,5],"Diagnostics",["BOUNDARY_CONDITIONS_CREATE_FINISH"])

        numberOfLoadIncrements = 1
        numberGlobalXElements = self.cantilever_elements[0]
        numberGlobalYElements = self.cantilever_elements[1]
        numberGlobalZElements = self.cantilever_elements[2]
        InterpolationType = 1
        if numberGlobalZElements == 0:
            numberOfXi = 2
        else:
            numberOfXi = 3

        # Get the number of computational nodes and this computational node number
        numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
        computationalNodeNumber = iron.ComputationalNodeNumberGet()

        # Create a 3D rectangular cartesian coordinate system
        self.coordinate_system = iron.CoordinateSystem()
        self.coordinate_system.CreateStart(coordinateSystemUserNumber)
        self.coordinate_system.DimensionSet(3)
        self.coordinate_system.CreateFinish()

        # Create a region and assign the coordinate system to the region
        self.region = iron.Region()
        self.region.CreateStart(regionUserNumber,iron.WorldRegion)
        self.region.LabelSet("Region 1")
        self.region.coordinateSystem = self.coordinate_system
        self.region.CreateFinish()

        # Define basis
        self.basis = iron.Basis()
        self.basis.CreateStart(basisUserNumber)
        if InterpolationType in (1,2,3,4):
            self.basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        elif InterpolationType in (7,8,9):
            self.basis.type = iron.BasisTypes.SIMPLEX
        self.basis.numberOfXi = numberOfXi
        self.basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfXi
        if NumberOfGaussXi>0:
            self.basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
        self.basis.CreateFinish()

        if UsePressureBasis:
            # Define pressure basis
            self.pressureBasis = iron.Basis()
            self.pressureBasis.CreateStart(pressureBasisUserNumber)
            if InterpolationType in (1,2,3,4):
                self.pressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
            elif InterpolationType in (7,8,9):
                self.pressureBasis.type = iron.BasisTypes.SIMPLEX
            self.pressureBasis.numberOfXi = numberOfXi
            self.pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_LAGRANGE]*numberOfXi
            if NumberOfGaussXi>0:
                self.pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
            self.pressureBasis.CreateFinish()

        # Start the creation of a generated mesh in the region
        self.generatedMesh = iron.GeneratedMesh()
        self.generatedMesh.CreateStart(generatedMeshUserNumber,self.region)
        self.generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
        if UsePressureBasis:
            self.generatedMesh.basis = [self.basis,self.pressureBasis]
        else:
            self.generatedMesh.basis = [self.basis]
        if numberGlobalZElements==0:
            self.generatedMesh.extent = [width,height]
            self.generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements]
        else:
            self.generatedMesh.extent = [width,length,height]
            self.generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]
        # Finish the creation of a generated mesh in the region
        self.mesh = iron.Mesh()
        self.generatedMesh.CreateFinish(meshUserNumber,self.mesh)

        # Create a decomposition for the mesh
        self.decomposition = iron.Decomposition()
        self.decomposition.CreateStart(decompositionUserNumber,self.mesh)
        self.decomposition.type = iron.DecompositionTypes.CALCULATED
        self.decomposition.CalculateFacesSet(True)
        self.decomposition.numberOfDomains = numberOfComputationalNodes
        self.decomposition.CreateFinish()

        # Create a field for the geometry
        self.geometricField = iron.Field()
        self.geometricField.CreateStart(geometricFieldUserNumber,self.region)
        self.geometricField.MeshDecompositionSet(self.decomposition)
        self.geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
        self.geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
        self.geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
        if InterpolationType == 4:
            self.geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        self.geometricField.CreateFinish()

        # Update the geometric field parameters from generated mesh
        self.generatedMesh.GeometricParametersCalculate(self.geometricField)

        # Create a fibre field and attach it to the geometric field
        self.fibreField = iron.Field()
        self.fibreField.CreateStart(fibreFieldUserNumber,self.region)
        self.fibreField.TypeSet(iron.FieldTypes.FIBRE)
        self.fibreField.MeshDecompositionSet(self.decomposition)
        self.fibreField.GeometricFieldSet(self.geometricField)
        self.fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
        if InterpolationType == 4:
            self.fibreField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        self.fibreField.CreateFinish()

        # Create the equations_set
        self.equationsSetField = iron.Field()
        self.equationsSet = iron.EquationsSet()
        self.equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
            iron.EquationsSetTypes.FINITE_ELASTICITY,
            iron.EquationsSetSubtypes.MOONEY_RIVLIN]
        self.equationsSet.CreateStart(equationsSetUserNumber,self.region,self.fibreField,
            self.equationsSetSpecification, equationsSetFieldUserNumber, self.equationsSetField)
        self.equationsSet.CreateFinish()

        # Create the dependent field
        self.dependentField = iron.Field()
        self.equationsSet.DependentCreateStart(dependentFieldUserNumber,self.dependentField)
        self.dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
        self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
        self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
        if UsePressureBasis:
            # Set the pressure to be nodally based and use the second mesh component
            if InterpolationType == 4:
                self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
                self.dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
            self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
            self.dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
        if InterpolationType == 4:
            self.dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        self.equationsSet.DependentCreateFinish()

        # Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
        iron.Field.ParametersToFieldParametersComponentCopy(
            self.geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
        iron.Field.ComponentValuesInitialiseDP(
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0)

        # Create the material fields.
        self.materialField = iron.Field()
        self.equationsSet.MaterialsCreateStart(materialFieldUserNumber,self.materialField)
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.V,"Density")

        # If the gel has two layers, use these lines to create the two component sets.
        self.materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U,1,iron.FieldInterpolationTypes.ELEMENT_BASED)
        self.materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U,2,iron.FieldInterpolationTypes.ELEMENT_BASED)

        self.equationsSet.MaterialsCreateFinish()

        # Set Mooney-Rivlin constants c10 and c01 respectively. Done in objective function using function.
        self.materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,density)

        # Create the source field with the gravity vector
        self.sourceField = iron.Field()
        self.equationsSet.SourceCreateStart(sourceFieldUserNumber,self.sourceField)
        if InterpolationType == 4:
            self.sourceField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
        else:
            self.sourceField.fieldScalingType = iron.FieldScalingTypes.UNIT
        self.equationsSet.SourceCreateFinish()

        # Set the gravity vector component values
        self.sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,gravity[0])
        self.sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,gravity[1])
        self.sourceField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,gravity[2])

        # Create equations
        self.equations = iron.Equations()
        self.equationsSet.EquationsCreateStart(self.equations)
        self.equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
        self.equations.outputType = iron.EquationsOutputTypes.NONE
        self.equationsSet.EquationsCreateFinish()

        # Define the problem
        self.problem = iron.Problem()
        self.problemSpecification = [iron.ProblemClasses.ELASTICITY,
                iron.ProblemTypes.FINITE_ELASTICITY,
                iron.ProblemSubtypes.NONE]
        self.problem.CreateStart(problemUserNumber, self.problemSpecification)
        self.problem.CreateFinish()

        # Create the problem control loop
        self.problem.ControlLoopCreateStart()
        self.controlLoop = iron.ControlLoop()
        self.problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],self.controlLoop)
        self.controlLoop.MaximumIterationsSet(numberOfLoadIncrements)
        self.problem.ControlLoopCreateFinish()

        # Create problem solver
        self.nonLinearSolver = iron.Solver()
        self.linearSolver = iron.Solver()
        self.problem.SolversCreateStart()
        self.problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,self.nonLinearSolver)

        if self.diagnostics == 4 or self.diagnostics == 'Matrix':
            self.nonLinearSolver.outputType = iron.SolverOutputTypes.MATRIX
        elif self.diagnostics == 3 or self.diagnostics == 'Solver':
            self.nonLinearSolver.outputType = iron.SolverOutputTypes.SOLVER
        elif self.diagnostics == 2 or self.diagnostics == 'Timing':
            self.nonLinearSolver.outputType = iron.SolverOutputTypes.TIMING
        elif self.diagnostics == 1 or self.diagnostics == 'Progress':
            self.nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
        else:
            self.nonLinearSolver.outputType = iron.SolverOutputTypes.NONE

        self.nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
        self.nonLinearSolver.NewtonAbsoluteToleranceSet(1e-10)
        self.nonLinearSolver.NewtonSolutionToleranceSet(1e-10)
        self.nonLinearSolver.NewtonRelativeToleranceSet(1e-10)
        self.nonLinearSolver.NewtonMaximumIterationsSet(int(1e6))
        self.nonLinearSolver.NewtonMaximumFunctionEvaluationsSet(int(1e6))
        self.nonLinearSolver.NewtonLinearSolverGet(self.linearSolver)
        self.linearSolver.linearType = iron.LinearSolverTypes.DIRECT
        #linearSolver.libraryType = iron.SolverLibraries.LAPACK
        self.problem.SolversCreateFinish()

        # Create solver equations and add equations set to solver equations
        self.solver = iron.Solver()
        self.solverEquations = iron.SolverEquations()
        self.problem.SolverEquationsCreateStart()
        self.problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,self.solver)
        self.solver.SolverEquationsGet(self.solverEquations)
        self.solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
        equationsSetIndex = self.solverEquations.EquationsSetAdd(self.equationsSet)
        self.problem.SolverEquationsCreateFinish()

        # Prescribe boundary conditions (absolute nodal parameters)
        self.boundaryConditions = iron.BoundaryConditions()
        self.solverEquations.BoundaryConditionsCreateStart(self.boundaryConditions)

        #for i in range(1, ((numberGlobalXElements+1)*(numberGlobalYElements+1)*(numberGlobalZElements+1)), (numberGlobalXElements+1)):
        #    self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        #    self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
        #    self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

        numberOfNodes = (NumberOfGaussXi + (NumberOfGaussXi-1)*(numberGlobalXElements-1))\
                        * (NumberOfGaussXi + (NumberOfGaussXi-1)*(numberGlobalYElements-1))\
                          * (NumberOfGaussXi + (NumberOfGaussXi-1)*(numberGlobalZElements-1))

        #for nodeNum in range(1, numberOfNodes):
        #    xValue = np.array([self.materialField.ParameterSetGetNodeDP(
        #        iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNum,1)])

            #if #np.isclose(xValue, 0.0):
        for nodeNum in range(1, numberOfNodes+1, (NumberOfGaussXi + (NumberOfGaussXi-1)*(numberGlobalXElements-1))):
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNum, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

        self.solverEquations.BoundaryConditionsCreateFinish()

    def set_Neo_Hookean_single_layer_parameter(self, parameter_value):
        """
        Call to change the material parameter values without executing all of the other calls in the setup function.

        :param parameter_value: The value of the c01 Neo-Hookean parameter.
        """

        numElements = self.cantilever_elements[0] * self.cantilever_elements[1] * self.cantilever_elements[2]

        for elementNum in range(1, numElements + 1):
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum, 1, parameter_value[0])
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum, 2, 0.0)

        iron.Field.ComponentValuesInitialiseDP(
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0.0)

    def set_Neo_Hookean_two_layer_parameters(self, parameter_values):
        """
        Call to set up a two layered cantilever after the rest of the cantilever has been set up.

        :param parameter_values: The value of the c01 Neo-Hookean parameter for each of the layers.
        """

        numElements = self.cantilever_elements[0] * self.cantilever_elements[1] * self.cantilever_elements[2]

        for elementNum in range(1, (numElements / 2) + 1):
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum, 1, parameter_values[0])
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum, 2, 0.0)
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum + (numElements / 2), 1, parameter_values[1])
            self.materialField.ParameterSetUpdateElementDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, elementNum + (numElements / 2), 2, 0.0)

        iron.Field.ComponentValuesInitialiseDP(
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0.0)

    def solve_simulation(self):
        self.problem.Solve()

    def export_results(self):
        self.fields = iron.Fields()
        self.fields.CreateRegion(self.region)
        self.fields.NodesExport("Cantilever","FORTRAN")
        self.fields.ElementsExport("Cantilever","FORTRAN")
        self.fields.Finalise()

    def generate_data(self, scale):
        """
        Create some artificial data which can be used to test if the optimisation routine can be used to find the
        parameters which generated the data.

        :param scale: Sets the number of points which will be generated.
                scale = 0 => Testing corner location mode
                scale = 1 => Determined by self.numPointsPerFace

        :return: Sets the simulation's data to a set of points on the surface of the FE model.
        """

        if scale == 0:
            dataLocations = np.array([[0, 0, 0]])

            elements = np.array([self.cantilever_elements[0], self.cantilever_elements[0] * self.cantilever_elements[1],
                                 self.cantilever_elements[0] + (self.cantilever_elements[0] * self.cantilever_elements[1]
                                 * (self.cantilever_elements[2] - 1)), self.cantilever_elements[0] * self.cantilever_elements[1] * self.cantilever_elements[2]])
            Xi = np.array([[1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1]])

            for i in range(4):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,elements[i],Xi[i],4)
                dataLocations = np.append(dataLocations, np.array([point[0:3]]), axis=0)

            dataLocations = dataLocations[1:, :]

        elif scale == 1:
            # First, select the relevant elements for each face.
            [leftEls, bottomEls, endEls, rightEls, topEls] = self.select_elements()

            # Next, create the grids of Xi values.
            [leftXi, bottomXi, endXi, rightXi, topXi] = self.create_Xi_grids()

            # Now use interpolation to find the coordinates at each Xi point.
            dataLocations = self.interpolate_Xi_points(leftEls, bottomEls, endEls, rightEls, topEls, leftXi, bottomXi, endXi, rightXi, topXi)

        return dataLocations

    def select_elements(self):
        """
        Finds which elements are on which face of the model and stores these. They can't be stored in a 2D array where
        each row corresponds to a particular face because each face can have a different number of elements on it, which
        does not allow for the formation of a 2D array.

        If this was really necessary, then the maximum number of elements could be found and a 2D array initialised with
        zeros on every row for that length. Then, for each row the elements could be copied into this array over the
        zeros. Since zero isn't a valid elements number, a loop could be used to extract the element numbers for each
        face from the corresponding row, stopping when the first zero is reached. This seems of limited use though.

        :return: Returns the sets of elements which have been calculated for each of the FE model's faces.
        """

        # Redefine the number of elements in each dimension for quicker use.
        x = self.cantilever_elements[0]
        y = self.cantilever_elements[1]
        z = self.cantilever_elements[2]

        # Select the correct elements for each face (put in own function)
        leftEls = rightEls = np.array([])
        for i in range(z):
            leftEls = np.append(leftEls, range(i*x*y+1, i*x*y+x+1))
            rightEls = np.append(rightEls, range(i*x*y+x*(y-1)+1, i*x*y+x*(y-1)+x+1))

        leftEls = leftEls.astype('int64')
        rightEls = rightEls.astype('int64')

        bottomEls = np.array(range(1, x*y+1))
        topEls = np.array(range((z-1)*x*y+1, x*y*z+1))
        endEls = np.array(range(x, x*y*z+1, x))

        return leftEls, bottomEls, endEls, rightEls, topEls

    def create_Xi_grids(self):
        """
        Using the specified number of points to be generated per side of each element face, the grids of Xi coordinates
        are created so that those points can be interpolated.

        :return: Returns the arrays of Xi coordinates for each of FE model's faces.
        """

        # With number of points per face defined, generate the Xi grids.
        array = np.linspace(0, 1, self.numPointsPerFace)
        [grid1, grid2] = np.meshgrid(array, array)

        for i in range(self.numPointsPerFace):
            for j in range(self.numPointsPerFace):
                if i == 0 and j == 0:
                    leftXi = np.array([[grid1[i,j], 0, grid2[i,j]]])
                    bottomXi = np.array([[grid1[i,j], grid2[i,j], 0]])
                    endXi = np.array([[1, grid1[i,j], grid2[i,j]]])
                    rightXi = np.array([[grid1[i,j], 1, grid2[i,j]]])
                    topXi = np.array([[grid1[i,j], grid2[i,j], 1]])

                else:
                    leftXi = np.append(leftXi, np.array([[grid1[i,j], 0, grid2[i,j]]]), axis=0)
                    bottomXi = np.append(bottomXi, np.array([[grid1[i,j], grid2[i,j], 0]]), axis=0)
                    endXi = np.append(endXi, np.array([[1, grid1[i,j], grid2[i,j]]]), axis=0)
                    rightXi = np.append(rightXi, np.array([[grid1[i,j], 1, grid2[i,j]]]), axis=0)
                    topXi = np.append(topXi, np.array([[grid1[i,j], grid2[i,j], 1]]), axis=0)

        return leftXi, bottomXi, endXi, rightXi, topXi

    def interpolate_Xi_points(self, leftEls, bottomEls, endEls, rightEls, topEls, leftXi, bottomXi, endXi, rightXi, topXi):
        """
        Using the specified elements for each face of the FE model and the set of Xi coordinates for that face, this
        function will interpolate the position of all the data points on each face of the FE model.

        :param leftEls: The set of elements for interpolating points on the model's left face.
        :param bottomEls: The set of elements for interpolating points on the model's bottom face.
        :param endEls: The set of elements for interpolating points on the model's end face.
        :param rightEls: The set of elements for interpolating points on the model's right face.
        :param topEls: The set of elements for interpolating points on the model's top face.
        :param leftXi: The set of Xi coordinates which specify where to interpolate points on the left face.
        :param bottomXi: The set of Xi coordinates which specify where to interpolate points on the bottom face.
        :param endXi: The set of Xi coordinates which specify where to interpolate points on the end face.
        :param rightXi: The set of Xi coordinates which specify where to interpolate points on the right face.
        :param topXi: The set of Xi coordinates which specify where to interpolate points on the top face.
        :return: Returns an N-by-3 array of data points at the specified coordinates. Note that there are no duplicates.
        """

        dataLocations = np.array([iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,leftXi[0],4)])

        for i in range(5):
            if i == 0:
                Xi = leftXi
                elements = leftEls
            elif i == 1:
                Xi = bottomXi
                elements = bottomEls
            elif i == 2:
                Xi = endXi
                elements = endEls
            elif i == 3:
                Xi = rightXi
                elements = rightEls
            elif i == 4:
                 Xi = topXi
                 elements = topEls

            for j in range(len(elements)):
                for k in range(len(Xi)):
                    point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,elements[j],Xi[k],4)
                    dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        from np_1_13_unique import unique
        dataLocations, dataIdx = unique(dataLocations[:,0:3], return_index=True, axis=0)
        dataLocations = dataLocations[np.argsort(dataIdx)]

        return dataLocations

    def point_projection(self):
        """
        Objective function for measuring the error between a FE model (contained in simulation) and the data measured from
        the surface of a real experiment.

        :return: error: A measure of the error between the data points and the surface of the FE model, stored as a float
                        in the simulation class itself.
        """

        errorValues = np.array([])
        errorVects = np.array([[0, 0, 0]]) # This must be initialised with a 1-by-3 array in order to allow appending of more such arrays
        dataIdx = 0

        # Projections must be done face by face. First step is to loop through the five faces onto which the data points
        # will be projected. Once the faces for projection are set for one group of elements, select the corresponding
        # points for projection and calculation their error before destroying both those projection and data sets to
        # make room for the next set.
        for i in range(2,7):
            numDataPoints, points, dataIdx = self.select_points(i, dataIdx)
            self.create_data_projection(i, numDataPoints, points)
            self.set_projection_elements(i)

            # Now perform the projections.
            projectionErr, projectionErrVec = self.perform_projection(numDataPoints)

            errorValues = np.append(errorValues, projectionErr)
            errorVects = np.append(errorVects, projectionErrVec, axis=0)

            # The data points stored at first by this function have now been passed to the data projection. To prepare
            # for the next loop of this calculation, now destroy the data points structure.
            self.dataPoints.Destroy()

        self.projection_error(errorValues, errorVects[1:,:]) # This removes the initialised 1-by-3 array so the data is not out of sync.

    def select_points(self, faceNum, prevDataIdx):
        """
        Finds which points to use for data projection depending on which face of the FE model is being used.

        :param faceNum: The OpenCMISS face ID to which the points will be projected, as an integer.
        :param prevDataIdx: The number of data points which have already been used for projections onto other faces,
                            as an integer
        :return:
        """

        # Set the number of elements in each direction to easier to use references.
        x = self.cantilever_elements[0]
        y = self.cantilever_elements[1]
        z = self.cantilever_elements[2]

        # To make the following algorithms actually readable, the number of points along a single edge is reassigned
        # to a smaller variable name.
        pts = self.numPointsPerFace

        # Now, depending on which face the projection_calculation function is up to, set the number of data points to
        # the number of elements on that face, multiplied by the number of points per element. Then also select that
        # sequence of data points out of the total data set.
        if faceNum == 2:
            numDataPoints = pts**2 + ((pts**2 - pts) * (z - 1)) + (((pts**2 - pts) + (z - 1) * (pts**2 - 2 * pts + 1))* (x - 1))
        elif faceNum == 3:
            numDataPoints = y * ((pts**2 - pts) + (pts**2 - 2 * pts + 1) * (x - 1))
        elif faceNum == 4:
            numDataPoints = y * z * (pts**2 - 2 * pts + 1)
        elif faceNum == 5:
            numDataPoints = x * z * (pts**2 - 2 * pts + 1)
        elif faceNum == 6:
            numDataPoints = x * ((y - 1) * (pts**2 - 2 * pts + 1) + (pts**2 - 3 * pts + 2))

        # Now isolate the points relevant to this face using the number of points calculated to exist on this face and
        # the record of how many points have already been used. Update how many points have been used as well.
        points = self.data[prevDataIdx:(prevDataIdx + numDataPoints)]
        prevDataIdx += numDataPoints

        return numDataPoints, points, prevDataIdx

    def create_data_projection(self, faceNum, numPoints, dataPoints):
        """
        Using OpenCMISS.iron, creates the necessary data projection class and applies all the relevant settings.

        :param faceNum: The OpenCMISS face ID to which the points will be projected, as an integer.
        :param numPoints: The number of points which are to be projected, as an integer.
        :param dataPoints: The coordinates of each data point to be projected, as an N-by-3 numpy array.
        """

        # The data points have been passed into this function, so loop through them after creating the structure to
        # house them so they can be used for projections.
        self.dataPoints = iron.DataPoints()
        self.dataPoints.CreateStart(self.region, numPoints)
        for pointNum, point in enumerate(dataPoints,1):
            self.dataPoints.ValuesSet(pointNum, point)
        self.dataPoints.CreateFinish()

        # Now create the data projection structure and pass in the data point structure which was just made.
        self.dataProjection = iron.DataProjection()
        self.dataProjection.CreateStart(faceNum, self.dataPoints, self.mesh)

        # Set tolerances and other settings for the data projection.
        self.dataProjection.AbsoluteToleranceSet(1.0e-15)
        self.dataProjection.RelativeToleranceSet(1.0e-15)
        self.dataProjection.MaximumNumberOfIterationsSet(int(1e9))
        self.dataProjection.ProjectionTypeSet(
            iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
        self.dataProjection.CreateFinish()

    def set_projection_elements(self, faceNum):
        """
        Selects the elements for the current projection, allowing the data points to be projected onto only the faces
        of the relevant elements of the FE model.

        :param faceNum: The OpenCMISS face ID of the face which will be projected onto, as an integer.
        """

        # Set the number of elements in each direction to easier to use references.
        x = self.cantilever_elements[0]
        y = self.cantilever_elements[1]
        z = self.cantilever_elements[2]

        # Set the elements for projection and the relevant faces. First check which iteration of the projcetion
        # projection calculation is current. Then select all the elements on the relevant side and set their outward
        # face to be the one to be used for projection.
        elements = np.array([])

        if faceNum == 2:
            for j in range(z):
                elements = np.append(elements, np.array(range((j*x*y + 1), (j*x*y + x + 1))))
        elif faceNum == 3:
            elements = np.array(range(1, (x*y + 1)))
        elif faceNum == 4:
            elements = np.array(range(x, x*y*z+1, x))
        elif faceNum == 5:
            for j in range(z):
                elements = np.append(elements, np.array(range((j*x*y + x*(y-1) + 1), ((j+1)*x*y + 1))))
        elif faceNum == 6:
            elements = np.array(range((x*y*(z-1) + 1), (x*y*z + 1)))

        elements = elements.astype('int32')
        faces = np.full(len(elements), faceNum, dtype='int32')
        self.dataProjection.ProjectionCandidatesSet(elements, faces)

    def perform_projection(self, numDataPoints):
        """
        Using the data projection structure previously set up in the simulation class, this function now calculates the
        distance (or the error) between the data points and the surface of the FE model which has most recently been
        solved.

        :param numDataPoints: The number of points which are to be evaluated, as an integer.
        :return:
        """

        self.dataProjection.DataPointsProjectionEvaluate(self.dependentField)
        projectionErr = np.zeros(numDataPoints)
        projectionErrVec = np.zeros((numDataPoints,3))
        for pointNum, pointIdx in enumerate(range(numDataPoints), 1):
            projectionErr[pointIdx] = self.dataProjection.ResultDistanceGet(pointNum)
            projectionErrVec[pointIdx,:] = self.dataProjection.ResultProjectionVectorGet(pointNum, 3)

        return projectionErr, projectionErrVec

    def projection_error(self, errorValues, errorVects):
        """
        Deals with the errors calculated from the data projections, both finding the RMS error between the two data sets
        as well as exporting the error vectors.

        :param errorValues: A 1D array containing the error values for each point in the two data sets.
        :param errorVects: A N-by-3 array containing the coordinate vectors which describe the difference between each
                            pair of points in the two data sets.
        """

        # Having found the projection errors, calculate the RMS error.
        for i in range(len(errorValues)):
            errorValues[i] = errorValues[i] ** 2
        self.error = np.sqrt(np.average(errorValues))

        # Now export the resultant vectors.
        exportDatapointsErrorExdata(self.data, errorVects, 'error', './', 'error')


def single_layer_objective_function(material_parameter, simulation):
    """
    The objective function is used in the optimal design and consists of setting the material parameters which were
    passed into the function, solving the FE model with those parameters and then projecting the data set onto the
    surface of the model to find the error between the two sets of material parameters.

    This can be used to determine what the true material parameters of an object are if the data set is obtained from
    an experiment with an unknown parameter. If the error output from this function is minimised, then the material
    parameter value which is being passed into this function must be close to that of the true parameter.

    :param material_parameter: The value of the material parameter, as a 1D numpy array with a float.
    :param simulation: A tuple containing the set up simulation.
    """

    simulation.set_Neo_Hookean_single_layer_parameter(material_parameter)
    simulation.solve_simulation()
    simulation.export_results()
    simulation.point_projection()

    return simulation.error

def two_layer_objective_function(material_parameters, simulation):
    """
    Another objective function used for parameter optimisation.

    :param material_parameters: The values of the material parameters, as a 1D numpy array of two floats.
    :param simulation: A tuple containing the set up simulation.

    """

    simulation.set_Neo_Hookean_two_layer_parameters(material_parameters)
    simulation.solve_simulation()
    simulation.export_results()
    simulation.point_projection()

    return simulation.error

###########
# Testing #
###########

if __name__ == "__main__":
    # Testing the use of the objective function.
    cantilever_dimensions = np.array([30, 12, 12])

    cantilever_elements = np.array([1, 1, 1])
    cantilever_true_parameter = np.array([1.2])
    cantilever_guess_parameter = np.array([1.2])

    cantilever_sim = CantileverSimulation()

    cantilever_sim.set_Xi_points_num(3)
    cantilever_sim.set_cantilever_dimensions(cantilever_dimensions)
    cantilever_sim.set_cantilever_elements(cantilever_elements)
    cantilever_sim.set_gravity_vector(np.array([0.0, 10, 0.0]))
    cantilever_sim.set_diagnostic_level(1)
    cantilever_sim.setup_cantilever_simulation()
    cantilever_sim.set_Neo_Hookean_single_layer_parameter(cantilever_true_parameter)
    cantilever_sim.solve_simulation()

    #data = cantilever_sim.generate_data(1)
    #print '1st Data Set'
    #print data
    #print '\n'
    #cantilever_sim.set_projection_data(data)

    cantilever_sim.set_projection_data()

    single_layer_objective_function(cantilever_guess_parameter, cantilever_sim)
    data2 = cantilever_sim.generate_data(1)[:,0:3]
    print '2nd Data Set'
    print '\n'
    print data2
    print '\n'
    print "RMS Error = {0}".format(cantilever_sim.error)
    print '\n'
