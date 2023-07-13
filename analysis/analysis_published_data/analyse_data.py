import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt

these_files = os.listdir(os.path.join(os.path.dirname(__file__),'csv_files'))
these_files.sort()
print(these_files)

column_data = []
for file in these_files:
    this_file = pd.read_csv(os.path.join(os.path.dirname(__file__),'csv_files',file))
    these_x_coordinates = this_file['X'].to_numpy()
    these_y_coordinates = this_file['Y'].to_numpy()
    these_positions = np.array([these_x_coordinates, these_y_coordinates]).transpose()
    column_data.append(these_positions)

rotated_columns = []
for index, column in enumerate(column_data):
    x_values = column[:,0]
    y_values = column[:,1]
    x_values -= x_values[0]
    y_values -= y_values[0]
    slope, offset = np.polyfit(x_values,y_values,1)
    reset_y_values = y_values - offset
    this_point_set = np.array((x_values,reset_y_values))
    angle = np.arctan(slope)
    rotation_angle = np.pi/2 - angle
    rotation_matrix = np.array([[np.cos(rotation_angle),-np.sin(rotation_angle)],
                                [np.sin(rotation_angle), np.cos(rotation_angle)]])
    rotated_points = rotation_matrix.dot(this_point_set)
    rotated_coordinates = rotated_points.transpose()
    rotated_x_values = rotated_coordinates[:,0]
    rotated_y_values = rotated_coordinates[:,1]
    scaling_factor = (len(x_values) - 1)/np.max(np.abs(rotated_y_values))
    x_values*=scaling_factor
    reset_y_values*=scaling_factor
    rotated_coordinates*=scaling_factor
    rotated_columns.append(rotated_coordinates)
    plt.figure()
    plt.scatter(x_values, reset_y_values, label = 'original')
    plt.scatter(rotated_x_values, rotated_y_values, label = 'rotated')
    plt.plot(x_values, slope*x_values)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.legend()
    plt.savefig(os.path.join(os.path.dirname(__file__),'output','column_' + str(index+1) + '_visualised.jpeg'))

envelope_projection_areas = []
column_lengths = []
for column in rotated_columns:
    x_values = column[:,0]
    difference = np.max(x_values) - np.min(x_values)
    column_length = len(x_values)
    column_lengths.append(column_length)
    envelope_projection_areas.append(difference**2)
    
envelope_projection_areas = np.array(envelope_projection_areas)
column_lengths = np.array(column_lengths)
print(envelope_projection_areas)
print('the average projection area for small columns')
short_envelope_projection_areas = envelope_projection_areas[column_lengths<6]
ea_mean = np.mean(short_envelope_projection_areas)
ea_std = np.std(short_envelope_projection_areas)
print('the mean is')
print(ea_mean)
print('the std is')
print(ea_std)
print('the number of cells is')
print(len(envelope_projection_areas))

plt.figure()
plt.scatter(column_lengths, envelope_projection_areas)
plt.xlabel('column length')
plt.ylabel('envelope projection area')
plt.savefig(os.path.join(os.path.dirname(__file__),'output','envelope_projection_areas.jpeg'))