from opencmiss.iron import iron
import numpy as np

def projectSurfaceData(self, data):
    """
    Project surface data onto model surface
    """
    numDatapoints = len(data)

    datapoints = iron.DataPoints()
    datapoints.CreateStart(self.simulation.region, numDatapoints)
    for pointNum, point in enumerate(data,1):
        datapoints.ValuesSet(pointNum, point)
    datapoints.CreateFinish()

    dataprojectionUserNumber = 1
    dataprojection = iron.DataProjection()
    dataprojection.CreateStart(dataprojectionUserNumber, datapoints,
                               self.simulation.mesh)
    dataprojection.AbsoluteToleranceSet(1.0e-15)
    dataprojection.RelativeToleranceSet(1.0e-15)
    dataprojection.MaximumNumberOfIterationsSet(int(1e9))
    dataprojection.ProjectionTypeSet(
        iron.DataProjectionProjectionTypes.BOUNDARY_FACES)
    dataprojection.CreateFinish()

    # \todo replace hardcoded elements with those on the top
    # surface of the gel.
    #            dataprojection.ProjectionCandidatesSet(
    #                range(1,25+1),[6]*25)

    dataprojection.DataPointsProjectionEvaluate(self.simulation.deformedField)
    # dataprojection.ProjectionEvaluate(equation_set.geometricField)
    # dataprojection.ProjectionEvaluate(equation_set.geometricField)
    projectionErr = np.zeros(numDatapoints)
    for pointNum, pointIdx in enumerate(range(numDatapoints), 1):
        projectionErr[pointIdx] = dataprojection.ResultDistanceGet(pointNum)
