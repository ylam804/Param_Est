from opencmiss.iron import iron
import numpy as np


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
        self.error = None
        self.fields = None

    def set_projection_data(self, data):
        """
        Initialises the set of data points which were scanned on the surface of the deformed experiment.

        :param data: A 2D array where each row is three point in mm of the x, y, and z coordinates respectively.
        """

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
        NumberOfGaussXi = 2

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
        # iron.DiagnosticsSetOn(iron.DiagnosticTypes.All,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

        numberOfLoadIncrements = 4
        numberGlobalXElements = self.cantilever_elements[0]
        numberGlobalYElements = self.cantilever_elements[1]
        numberGlobalZElements = self.cantilever_elements[2]
        InterpolationType = 1
        if numberGlobalZElements==0:
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
        self.basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
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
            self.pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
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
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,-8.0)

        # Create the material fields.
        self.materialField = iron.Field()
        self.equationsSet.MaterialsCreateStart(materialFieldUserNumber,self.materialField)
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.V,"Density")
        self.equationsSet.MaterialsCreateFinish()

        ##############################################################
        ##                                                          ##
        ## Need to allow for different constitutive relation and/or ##
        ##          parameter values to be selected.                ##
        ##                                                          ##
        ##############################################################

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

        ##############################################################
        ##                                                          ##
        ## Need to automate BC prescribing for any sized cantilever.##
        ##                                                          ##
        ##############################################################

        # Prescribe boundary conditions (absolute nodal parameters)
        self.boundaryConditions = iron.BoundaryConditions()
        self.solverEquations.BoundaryConditionsCreateStart(self.boundaryConditions)

        for i in range(1, (numberGlobalXElements*numberGlobalYElements*numberGlobalZElements), (numberGlobalXElements+1)):
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

        # Set x=0 nodes to no x displacement
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,1,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,3,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,5,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,7,1,iron.BoundaryConditionsTypes.FIXED,0.0)

        # Set y=0 nodes to no y displacement
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,1,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,3,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,5,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,7,2,iron.BoundaryConditionsTypes.FIXED,0.0)

        # Set z=0 nodes to no y displacement
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,1,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,3,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,5,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        #self.boundaryConditions.AddNode(self.dependentField,iron.FieldVariableTypes.U,1,1,7,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        self.solverEquations.BoundaryConditionsCreateFinish()

    def set_Mooney_Rivlin_parameter_values(self, parameter_values):
        """
        Call to change the material parameter values without executing all of the other calls in the setup function.

        :param parameter_values: The values of the c01 and c10 Mooney-Rivlin parameters respectively.
        """

        self.materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1, parameter_values[0])
        self.materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1, parameter_values[1])

    def solve_simulation(self):
        self.problem.Solve()

    def export_results(self):
        self.fields = iron.Fields()
        self.fields.CreateRegion(self.region)
        self.fields.NodesExport("Cantilever","FORTRAN")
        self.fields.ElementsExport("Cantilever","FORTRAN")
        self.fields.Finalise()

    def prepare_projection(self):
        """
        Prepare the simulation object for data projection processes. This should only be run once as the dataPoints and
        dataProjection structures only need to be made once.

        :param simulation: A set up FE model.
        :param data: An array of length n where each row is 3D the coordinates of a single data point measured from the
                        surface of the real experiment.
        """

        numDatapoints = len(self.data)

        self.dataPoints = iron.DataPoints()
        self.dataPoints.CreateStart(self.region, numDatapoints)
        for pointNum, point in enumerate(self.data,1):
            self.dataPoints.ValuesSet(pointNum, point)
        self.dataPoints.CreateFinish()

        dataProjectionUserNumber = 1
        self.dataProjection = iron.DataProjection()
        self.dataProjection.CreateStart(dataProjectionUserNumber, self.dataPoints,
                                        self.mesh)
        # Set tolerances and other settings for the data projection.
        self.dataProjection.AbsoluteToleranceSet(1.0e-15)
        self.dataProjection.RelativeToleranceSet(1.0e-15)
        self.dataProjection.MaximumNumberOfIterationsSet(int(1e9))
        self.dataProjection.ProjectionTypeSet(
            iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
        self.dataProjection.CreateFinish()

        # Set the elements for projection and the relevant faces.
        pos = 1
        elements = np.array([])
        faces = np.array([])
        for z in range(self.cantilever_elements[2]):
            for y in range(self.cantilever_elements[1]):
                for x in range(self.cantilever_elements[0]):
                    # Check if the element should be a candidate for projection.
                    if x == 0 or x == (self.cantilever_elements[0]-1) or y == 0 or y == (self.cantilever_elements[1]-1) or z == 0 or z == (self.cantilever_elements[2]-1):
                        elements = np.append(elements, np.array([pos]))
                        faces = np.append(faces, np.array([4]))
                        # Now find which of its faces should be projected onto.

                        ## Info on face locations:
                        #
                        # Coordinates at (0,0,0), no error on face 1
                        # Coordinates at (0,20,20), no error => face 1 is base x-direction
                        #
                        # Coordinates at (0,0,0), no error on face 2
                        # Coordinates at (30,0,20), error = 0.02 => face 2 is base y-direction
                        #
                        # Coordinates at (30,20,-3), error = 0.10 on face 3
                        # Coordinates at (30,20,37), error = 39.9 => face 3 is base z-direction
                        #
                        # Coordinates at (60,0,0), error = 2.04 on face 4
                        # Coordinates at (58.5,0,0), error = 0.5 on face 4
                        # Coordinates at (58.5,20,20), error = 2.12 => face 4 is end of x-direction
                        #
                        # Coordinates at (30,40,20), error = 0.02 on face 5
                        # Coordinates at (40,40,10), error = 0.08 => face 5 is end of y-direction
                        #
                        # Coordinates at (30,20,37), error = 0.03 on face 6
                        # Coordinates at (40,10,37), error = 0.95 => face 6 is end of z-direction


                        #if x == 0:
                        #    faces = np.append(faces, np.array([1]))
                        #elif x == (elementArray[0]-1):
                        #    faces = np.append(faces, np.array([6]))
                        #elif y == (elementArray[1]-1):
                        #    faces = np.append(faces, np.array([2]))
                        #elif y == (elementArray[1]-1):
                        #    faces = np.append(faces, np.array([5]))
                        #elif z == (elementArray[2]-1):
                        #    faces = np.append(faces, np.array([3]))
                        #else:
                        #    faces = np.append(faces, np.array([4]))
                    # Increase the position for the next iteration.
                    pos += 1

        elements = elements.astype('int32')
        faces = faces.astype('int32')
        self.dataProjection.ProjectionCandidatesSet(elements, faces)

    def projection_calculation(self):
        """
        Objective function for measuring the error between a FE model (contained in simulation) and the data measured from
        the surface of a real experiment.

        :param simulation: A set up FE model.
        :param data: An array of length n where each row is 3D the coordinates of a single data point measured from the
                        surface of the real experiment.
        :return: error: a measure of the error between the data points and the surface of the FE model.
        """

        numDatapoints = len(self.data)

        # Now perform the projections.
        self.dataProjection.DataPointsProjectionEvaluate(self.dependentField)
        projectionErr = np.zeros(numDatapoints)
        for pointNum, pointIdx in enumerate(range(numDatapoints), 1):
            projectionErr[pointIdx] = self.dataProjection.ResultDistanceGet(pointNum)

        # Having found the projection errors, calculate the RMS error.
        self.error = 0
        for i in range(numDatapoints):
            self.error += projectionErr[i] ** 2

        self.error = np.sqrt(np.average(self.error))

        return self.error

    #def generate_data(self):
    #    ParameterSetInterpolateSingleXiDP(self, variableType, fieldSetType, derivativeNumber, userElementNumber, xi, valuesSize)


def cantilever_objective_function(x, cantilever_simulation):
    cantilever_simulation.set_Mooney_Rivlin_parameter_values(x)
    cantilever_simulation.solve_simulation()
    cantilever_simulation.export_results()
    cantilever_simulation.error = cantilever_sim.projection_calculation()

    return cantilever_simulation.error

###########
# Testing #
###########


# Testing the use of the objective function.
data = np.array([[58, 0, 0], [58, 40, 0], [58, 0, 40], [58, 40, 40], [58, 20, 20]])
cantilever_dimensions = np.array([60, 40, 40])
cantilever_elements = np.array([1, 1, 1])
cantilever_initial_parameters = np.array([1, 1])

cantilever_sim = CantileverSimulation()
cantilever_sim.set_projection_data(data)
cantilever_sim.set_cantilever_dimensions(cantilever_dimensions)
cantilever_sim.set_cantilever_elements(cantilever_elements)
cantilever_sim.set_diagnostic_level(0)
cantilever_sim.setup_cantilever_simulation()
cantilever_sim.prepare_projection()

error = cantilever_objective_function(cantilever_initial_parameters, cantilever_sim)
print error
print '\n\n'

#cantilever_second_parameters = np.array([2, 1])
#error = cantilever_objective_function(cantilever_second_parameters, cantilever_sim)
#print error
#print '\n\n'

#cantilever_third_parameters = np.array([1, 2])
#error = cantilever_objective_function(cantilever_third_parameters, cantilever_sim)
#print error
