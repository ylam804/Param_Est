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

def exportDataointsHDF5(data, error, label, directory, filename):

    import h5py
    hdf5_filename = '{0}/{1}.h5'.format(directory, filename)
    hdf5_main_grip = h5py.File(hdf5_filename, 'w')

    hdf5_main_grip.create_dataset(label, data)

    hdf5_main_grip.close()

def exportDatapointsErrorExdata(data, error, label, directory,
                                   filename):
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
