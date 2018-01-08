import numpy as np

#
class ModelFace:
    """
    A class which represents one side of a FE simulation
    """

    def __init__(self):
        """
        Initialise a new face.
        """

        self.data = np.array([])
        self.OpenCMISSFaceID = None
        self.OpenCMISSElements = None

    def add_data(self, data):
        """
        Add a list of data points to this face.

        :param data: An N-by-3 numpy array of data points, with each row representing a single point as the x, y and z
                        coordinates in that order.
        """

        self.data = data

    def set_face_ID(self, faceID):
        """
        Set the face number which is used by OpenCMISS for projecting the data points onto.

        :param faceID: The OpenCMISS face ID as an integer.
        """

        self.OpenCMISSFaceID = faceID

    def set_elements(self, elements):
        """
        Set the element numbers which will be used by OpenCMISS for projection of the data points onto.

        :param elements: A numpy array of integers which represent the elements for projection using the given face.
        :return:
        """
