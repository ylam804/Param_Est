from opencmiss.iron import iron
import numpy as np
import math


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
        k = np.array([0, 0, 1])

        # Calculate the cross of k and v
        kcrossV = np.array([(k[1]*v[2] - k[2]*v[1]), (k[2]*v[0] - k[0]*v[2]), (k[0]*v[1] - k[1]*v[0])])

        # Now calculate the dot of k and v
        kdotV = k[0]*v[0] + k[1]*v[1] + k[2]*v[2]

        # Now in the final calculation we replace the angel theta in the formula with phi because that is the angle
        # which denotes the rotation of the beam about its longitudinal axis.
        v = v * math.cos(phi) + (kcrossV * math.sin(phi)) + (k * kdotV * (1 - math.cos(phi)))

        return v

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
            self.dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,0)

        # Create the material fields.
        self.materialField = iron.Field()
        self.equationsSet.MaterialsCreateStart(materialFieldUserNumber,self.materialField)
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
        self.materialField.VariableLabelSet(iron.FieldVariableTypes.V,"Density")


        #loop through component number
        #    self.materialField.ComponentInterpolationSet(iron.FieldVariableTypes.U,componentNumber,iron.FieldInterpolationTypes.ELEMENT_BASED)


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


        #self.materialField.ParameterSetUpdateElementDP(variableType, fieldSetType, userElementNumber, componentNumber, value)

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

        # Prescribe boundary conditions (absolute nodal parameters)
        self.boundaryConditions = iron.BoundaryConditions()
        self.solverEquations.BoundaryConditionsCreateStart(self.boundaryConditions)

        for i in range(1, ((numberGlobalXElements+1)*(numberGlobalYElements+1)*(numberGlobalZElements+1)), (numberGlobalXElements+1)):
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
            self.boundaryConditions.AddNode(self.dependentField, iron.FieldVariableTypes.U, 1, 1, i, 3, iron.BoundaryConditionsTypes.FIXED, 0.0)

        self.solverEquations.BoundaryConditionsCreateFinish()

    def set_Mooney_Rivlin_parameter_values(self, parameter_values):
        """
        Call to change the material parameter values without executing all of the other calls in the setup function.

        :param parameter_values: The values of the c01 and c10 Mooney-Rivlin parameters respectively.
        """

        self.materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1, parameter_values[0])
        self.materialField.ComponentValuesInitialiseDP(
            iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 2, 0.0)
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

    def prepare_projection(self, i):
        """
        Prepare the simulation object for data projection processes. This should only be run once as the dataPoints and
        dataProjection structures only need to be made once.

        :param simulation: A set up FE model.
        :param data: An array of length n where each row is 3D the coordinates of a single data point measured from the
                        surface of the real experiment.
        """

        # Set the number of elements in each direction to easier to use references.
        x = self.cantilever_elements[0]
        y = self.cantilever_elements[1]
        z = self.cantilever_elements[2]

        size = 11

        # Now, depending on which face the projection_calculation function is up to, set the number of data points to
        # the number of elements on that face, multiplied by the number of points per element. Then also select that
        # sequence of data points out of the total data set.
        if i == 2 or i == 5:
            numDataPoints = x * z * size**2
            if i == 2:
                points = self.data[1:numDataPoints+1, 0:3]
            elif i == 5:
                points = self.data[((x*z + x*y + y*z)*size**2 + 1):((2*x*z + x*y + y*z)*size**2 + 1), 0:3]
        elif i == 3 or i == 6:
            numDataPoints = x * y * size**2
            if i == 3:
                points = self.data[(x*z*size**2 + 1):(x*z + x*y)*size**2+ 1, 0:3]
            elif i == 6:
                points = self.data[((2*x*z + x*y + y*z)*size**2 + 1):((2*x*z + 2*x*y + y*z)*size**2 + 1), 0:3]
        elif i == 4:
            numDataPoints = y * z * size**2
            points = self.data[((x*z + x*y)*size**2 + 1):((x*z + x*y + y*z)*size**2 + 1), 0:3]

        # Having defined which points are to be used in this pass of the projection calculation, now create the data
        # point structure.
        self.dataPoints = iron.DataPoints()
        self.dataPoints.CreateStart(self.region, numDataPoints)
        for pointNum, point in enumerate(points,1):
            self.dataPoints.ValuesSet(pointNum, point)
        self.dataPoints.CreateFinish()

        # Now create the data projection structure and pass in the data point structure which was just made.
        dataProjectionUserNumber = i
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

        # Set the elements for projection and the relevant faces. First check which iteration of the projcetion
        # projection calculation is current. Then select all the elements on the relevant side and set their outward
        # face to be the one to be used for projection.
        elements = np.array([])

        if i == 2:
            for j in range(z):
                elements = np.append(elements, np.array(range((j*x*y + 1), (j*x*y + x + 1))))
        elif i == 3:
            elements = np.array(range(1, (x*y + 1)))
        elif i == 4:
            elements = np.array(range(x, x*y*z+1, x))
        elif i == 5:
            for j in range(z):
                elements = np.append(elements, np.array(range((j*x*y + x*(y-1) + 1), ((j+1)*x*y + 1))))
        elif i == 6:
            elements = np.array(range((x*y*(z-1) + 1), (x*y*z + 1)))

        elements = elements.astype('int32')
        faces = np.full(len(elements), i, dtype='int32')
        self.dataProjection.ProjectionCandidatesSet(elements, faces)
        return numDataPoints

    def projection_calculation(self):
        """
        Objective function for measuring the error between a FE model (contained in simulation) and the data measured from
        the surface of a real experiment.

        :return: error: a measure of the error between the data points and the surface of the FE model.
        """

        errorValues = np.array([])

        # Projections must be done face by face. First step is to loop through the five faces onto which the data points
        # will be projected. Once the faces for projection are set for one group of elements, select the corresponding
        # points for projection and calculation their error before destroying both those projection and data sets to
        # make room for the next set.
        for i in range(2,7):
            numDataPoints = self.prepare_projection(i)

            # Now perform the projections.
            self.dataProjection.DataPointsProjectionEvaluate(self.dependentField)
            projectionErr = np.zeros(numDataPoints)
            projectionErrVec = np.zeros((numDataPoints,3))
            for pointNum, pointIdx in enumerate(range(numDataPoints), 1):
                projectionErr[pointIdx] = self.dataProjection.ResultDistanceGet(pointNum)
                projectionErrVec[pointIdx,:] = self.dataProjection.ResultProjectionVectorGet( pointNum, 3)

            errorValues = np.append(errorValues, projectionErr)



            # The data points stored at first by this function have now been passed to the data projection. To prepare
            # for the next loop of this calculation, now destroy the data points structure.
            self.dataPoints.Destroy()

        exportDatapointsErrorExdata(self.data[726*3+484:,:], projectionErrVec, 'error', './', 'error')

        # Having found the projection errors, calculate the RMS error.
        self.error = 0
        for i in range(len(errorValues)):
            errorValues[i] = errorValues[i] ** 2

        self.error = np.sqrt(np.average(errorValues))

        return self.error

    def generate_data(self, scale):
        """
        Create some artificial data which can be used to test if the optimisation routine can be used to find the
        parameters which generated the data.

        :param scale: Sets the number of points which will be generated.
                Coarse = 0
                Moderate = 1
                Fine = 2
                Testing = 3     ** Note: model should only contain one element **
        :return: Sets the simulation's data to a set of points on the surface of the FE model.
        """

        # Redefine the number of elements in each dimension for quicker use.
        x = self.cantilever_elements[0]
        y = self.cantilever_elements[1]
        z = self.cantilever_elements[2]

        dataLocations = np.array([])

        leftElems = rightElems = np.array([])
        for i in range(z):
            leftElems = np.append(leftElems, range(i*x*y+1, i*x*y+x+1))
            rightElems = np.append(rightElems, range(i*x*y+x*(y-1)+1, i*x*y+x*(y-1)+x+1))

        leftElems = leftElems.astype('int64')
        rightElems = rightElems.astype('int64')

        bottomElems = np.array(range(1, x*y+1))
        topElems = np.array(range((z-1)*x*y+1, x*y*z+1))
        endElems = np.array(range(x, x*y*z+1, x))

        size = 11

        array = np.linspace(0, 1, size)
        [grid1, grid2] = np.meshgrid(array, array)

        for i in range(size):
            for j in range(size):
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

        dataLocations = np.array([iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,leftXi[0],4)])

        for i in range(len(leftElems)):
            for j in range(len(leftXi)):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,leftElems[i],leftXi[j],4)
                dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        for i in range(len(rightElems)):
            for j in range(len(bottomXi)):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,bottomElems[i],bottomXi[j],4)
                dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        for i in range(len(endElems)):
            for j in range(len(endXi)):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,endElems[i],endXi[j],4)
                dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        for i in range(len(rightElems)):
            for j in range(len(rightXi)):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,rightElems[i],rightXi[j],4)
                dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        for i in range(len(topElems)):
            for j in range(len(topXi)):
                point = iron.Field_ParameterSetInterpolateSingleXiDPNum(1,4,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,topElems[i],topXi[j],4)
                dataLocations = np.append(dataLocations, np.array([point]), axis=0)

        # Note: these loops could be reduced down to just three since the left and right as well as the top and bottom
        #       could be looped through at the same time. This hasn't been done so that all the points on a single face
        #       of the model are grouped consecutively together in the dataLocations array.

        return dataLocations

def cantilever_objective_function(x, simulation):

    simulation.set_Mooney_Rivlin_parameter_values(x)
    simulation.solve_simulation()
    simulation.export_results()
    simulation.error = simulation.projection_calculation()

    return simulation.error

def exportDatapointsErrorExdata(data, error, label, directory,
                                   filename):

    data = data[1:,:]
    # Shape of data should be a [num_datapoints,dim] numpy array.
    field_id = open(directory + filename + '.exdata', 'w')
    field_id.write(' Group name: {0}\n'.format(label))
    field_id.write(' #Fields=3\n')
    field_id.write(
        ' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    field_id.write('   x.  Value index= 1, #Derivatives=0\n')
    field_id.write('   y.  Value index= 2, #Derivatives=0\n')
    field_id.write('   z.  Value index= 3, #Derivatives=0\n')
    field_id.write(' 2) error, field, rectangular cartesian, #Components=3\n')
    field_id.write('   x.  Value index= 4, #Derivatives=0\n')
    field_id.write('   y.  Value index= 5, #Derivatives=0\n')
    field_id.write('   z.  Value index= 6, #Derivatives=0\n')
    field_id.write(' 3) scale, field, rectangular cartesian, #Components=3\n')
    field_id.write('   x.  Value index= 7, #Derivatives=0\n')
    field_id.write('   y.  Value index= 8, #Derivatives=0\n')
    field_id.write('   z.  Value index= 9, #Derivatives=0\n')

    for point_idx, point in enumerate(range(1, data.shape[0] + 1)):
        field_id.write(' Node: {0}\n'.format(point))
        for value_idx in range(data.shape[1]):
            field_id.write(' {0:.12E}\n'.format(data[point_idx, value_idx]))
        for value_idx in range(data.shape[1]):
            field_id.write(' {0:.12E}\n'.format(error[point_idx, value_idx]))
        # Scale field is absolute of the error field to ensure vector
        # direction does not change.
        for value_idx in range(data.shape[1]):
            field_id.write(
                ' {0:.12E}\n'.format(abs(error[point_idx, value_idx])))
    field_id.close()


def exportDatapointsExdata(data, label, directory, filename):
    # Shape of data should be a [num_datapoints,dim] numpy array.

    field_id = open(directory + filename + '.exdata', 'w')
    field_id.write(' Group name: {0}\n'.format(label))
    field_id.write(' #Fields=1\n')
    field_id.write(
        ' 1) coordinates, coordinate, rectangular cartesian, #Components=3\n')
    field_id.write('  x.  Value index=1, #Derivatives=0, #Versions=1\n')
    field_id.write('  y.  Value index=2, #Derivatives=0, #Versions=1\n')
    field_id.write('  z.  Value index=3, #Derivatives=0, #Versions=1\n')

    for point_idx, point in enumerate(range(1, data.shape[0] + 1)):
        field_id.write(' Node: {0}\n'.format(point))
        for value_idx in range(data.shape[1]):
            field_id.write(' {0:.12E}\n'.format(data[point_idx, value_idx]))
    field_id.close()

    import h5py
    hdf5_filename = '{0}/{1}.h5'.format(directory, filename)
    hdf5_main_grip = h5py.File(hdf5_filename, 'w')

    hdf5_main_grip.create_dataset(label, data)

    hdf5_main_grip.close()

###########
# Testing #
###########

if __name__ == "__main__":
    # Testing the use of the objective function.
    data = np.array([[54.127, 0.724, -11.211], [54.127, 39.276, -11.211], [64.432, -0.669, 27.737], [64.432, 40.669, 27.737]])
    cantilever_dimensions = np.array([30, 12, 12])
    cantilever_elements = np.array([3, 2, 2])
    cantilever_initial_parameter = np.array([2.05])
    cantilever_guess_parameter = np.array([2.05])

    cantilever_sim = CantileverSimulation()

    cantilever_sim.set_cantilever_dimensions(cantilever_dimensions)
    cantilever_sim.set_cantilever_elements(cantilever_elements)
    cantilever_sim.set_gravity_vector(np.array([0.0, 0.0, -9.81]))
    cantilever_sim.set_diagnostic_level(1)
    cantilever_sim.setup_cantilever_simulation()
    cantilever_sim.solve_simulation()

    data = cantilever_sim.generate_data(3)[:,0:3]

    #from input-output import exportDatapointsExdata
    #exportDatapointsExdata(data, 'test', '/hpc/jada201/opt/CantileverParameterOptimisation/', 'data')

    print data
    cantilever_sim.set_projection_data(data)

    error = cantilever_objective_function(cantilever_guess_parameter, cantilever_sim)
    print('RMS Error = ', error)
    print '\n'
