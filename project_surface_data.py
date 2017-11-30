from opencmiss.iron import iron
import numpy as np


def projectSurfaceData(self, data, elementArray):
    """
    Project surface data onto model surface
    """

    #Find the number of data points which need to be projected.
    numDatapoints = len(data)

    # Start making the iron.DataPoints object which is used in the projection.
    datapoints = iron.DataPoints()
    # Creates the DataPoints object, requiring a region (which I don't understand where it comes from, what is the
    # simulation bit?) and the number of data points.
    datapoints.CreateStart(self.simulation.region, numDatapoints)
    # Using enumerate, the values from the data are moved into the DataPoints object as the PointNum variable
    # counts which point the loop is up to while the point variable extracts the information out of the provided data.
    for pointNum, point in enumerate(data,1):
        datapoints.ValuesSet(pointNum, point)
    datapoints.CreateFinish()

    # Start setting up the projection.
    dataprojectionUserNumber = 1
    dataprojection = iron.DataProjection()
    # Starting the DataProjection object requires the DataPoints object which was just defined above, as well as a mesh
    # (which I also don't understand because it uses the simulation?).
    dataprojection.CreateStart(dataprojectionUserNumber, datapoints,
                               self.simulation.mesh)
    # Set tolerances and other settings for the data projection.
    dataprojection.AbsoluteToleranceSet(1.0e-15)
    dataprojection.RelativeToleranceSet(1.0e-15)
    dataprojection.MaximumNumberOfIterationsSet(int(1e9))
    dataprojection.ProjectionTypeSet(
        iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
    dataprojection.CreateFinish()

    # Now select the set of elements which will be possible candidates for the data projection. At the same time, chose
    # which faces of each element will be used to project onto.
    pos = 0
    elements = np.array([])
    faces = np.array([])
    for z in range(elementArray[2]):
        for y in range(elementArray[1]):
            for x in range(elementArray[0]):
                # Check if the element should be a candidate for projection.
                if x == 0 or x == (elementArray[0]-1) or y == 0 or y == (elementArray[1]-1) or z == 0 or z == (elementArray[2]-1):
                    elements = np.append(elements, np.array([pos]))
                    # Now find which of its faces should be projected onto.
                    if x == 0:
                        faces = np.append(faces, np.array([1]))
                    elif x == (elementArray[0]-1):
                        faces = np.append(faces, np.array([6]))
                    elif y == (elementArray[1]-1):
                        faces = np.append(faces, np.array([2]))
                    elif y == (elementArray[1]-1):
                        faces = np.append(faces, np.array([5]))
                    elif z == (elementArray[2]-1):
                        faces = np.append(faces, np.array([3]))
                    else:
                        faces = np.append(faces, np.array([4]))
                # Increase the position for the next iteration.
                pos += 1

    dataprojection.ProjectionCandidatesSet(elements, faces)

    # Now perform the projections.
    dataprojection.DataPointsProjectionEvaluate(self.simulation.deformedField)
    # dataprojection.ProjectionEvaluate(equation_set.geometricField)
    projectionErr = np.zeros(numDatapoints)
    for pointNum, pointIdx in enumerate(range(numDatapoints), 1): # Will this just overwrite the first index of projectionErr over and over?
        projectionErr[pointIdx] = dataprojection.ResultDistanceGet(pointNum)

    # Having found the projection errors, calculate the RMS error.
    RMSError = 0
    for i in range(numDatapoints):
        RMSError += projectionErr[i] ** 2

    RMSError = np.sqrt(np.average(RMSError))

    return RMSError


#def
#
#    nodal_values = np.zeros((len(list_of_nodes), 3))
#    # Set the geometric information from the exnode file
#    # Read the geometric field
#    for node_idx, node_num in enumerate(list_of_nodes):
#        version = 1
#        for component in range(1, numberOfXi + 1):
#            if interpolation == "linear" or interpolation == "cubicLagrange":
#                derivs = [1]
#           elif interpolation == "cubicHermite":
#               derivs = range(1, 9)
#           for derivative in derivs:
#                nodal_values[node_idx, component - 1] = geometricField.ParameterSetGetNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 1, derivative, node_num, component)

###########
# TESTING #
###########


elements = np.array([])
numXElements = 4
numYElements = 8
numZElements = 4
pos = 0

#for x in range(numXElements):
#    for y in range(numYElements):
#        for z in range(numZElements):
#            if x == 0 or x == (numXElements-1) or y == 0 or y == (numYElements-1) or z == 0 or z == (numZElements-1):
#                elements = np.append(elements, np.array([pos]))
#            pos += 1

#print elements

data = np.array([1, 2, 3])

error = projectSurfaceData(data, numXElements, numYElements, numZElements)
