# -*- coding: utf-8 -*-
"""
Useful methods for calculating metrics on the cartilage sheet simulations generated from Chaste. 

Created on Thu Aug 10 16:14:20 2017

@author: Sonja Mathias
"""
import numpy as np
import itertools
import matplotlib.pyplot as plt


def read_in_data(path_to_dir):
    """
    
    """
    content = [];
    with open(path_to_dir+'results.viznodes') as f:
        for line in f:
            split_line = line.split()
            content.append(split_line)

    data = np.array(list(itertools.zip_longest(*content, fillvalue=np.nan)), dtype=float).T
    return (data[:,0], data[:,1:])
    
def get_border_layers(path_to_dir):
    
    # get cell tissue types
    content_ctt = [];
    with open(path_to_dir + 'results.vizcelltissuetypes') as f:
        for line in f:
            split_line = line.split()
            content_ctt.append(split_line)
    data_ctt = np.array(list(itertools.zip_longest(*content_ctt, fillvalue=np.nan)), dtype=float).T 

    cell_tissue_types = data_ctt[:,1:]
    


    # get cell division directions
    content_cdd = [];
    with open(path_to_dir + 'results.vizcelldivisiondirections') as f:
        for line in f:
            split_line = line.split()
            content_cdd.append(split_line)
    data_cdd = np.array(list(itertools.zip_longest(*content_cdd, fillvalue=np.nan)), dtype=float).T   
    
    cell_division_directions = data_cdd[:,1:]
    
    #print(np.sum((cell_tissue_types == 23) & (cell_division_directions == 17) , axis = 1))
    
    # get mask for lower boundary layer
    mask_lower = (cell_tissue_types ==  23) & (cell_division_directions == 17)
    mask_lower = np.kron(mask_lower, np.array([True,True,True]))
    # get mask for upper boundary layer
    mask_upper = (cell_tissue_types ==  23) & (cell_division_directions == 18)
    mask_upper = np.kron(mask_upper, np.array([True,True,True]))
    
    return (mask_lower, mask_upper)
    
def get_clonal_patches(path_to_dir):
    
    #print('Extracting clonal patches')
    # get the ancestors
    content = [];
    with open(path_to_dir + 'results.vizancestors') as f:
        for line in f:
            split_line = line.split()
            content.append(split_line)
    data = np.array(list(itertools.zip_longest(*content, fillvalue=np.nan)), dtype=float).T 

    ancestors = data[:,1:]

    # get all possible ancestor ids from the first line, -1 means not set
    ids = np.unique(ancestors[0, ~np.isnan(ancestors[0,:])]) #only search over non-nan elements

    # make a boolean mask for each of them
    clonal_patches = {}
    for patch in ids:
        mask = ancestors == patch
        mask = np.kron(mask, np.array([True,True,True]))
        clonal_patches[patch] = mask  
    
    # return a map of ancestor indices and corresponding boolean masks
    return clonal_patches

def calculate_width(coordinates):

    n_times, n_nodes = coordinates.shape
    width =  np.zeros((n_times, 3));
    width[:, 0] = abs(np.nanmax(coordinates[:,0::3], axis=1)- np.nanmin(coordinates[:,0::3], axis=1))
    width[:, 1] = abs(np.nanmax(coordinates[:,1::3], axis=1)- np.nanmin(coordinates[:,1::3], axis=1))
    width[:, 2] = abs(np.nanmax(coordinates[:,2::3], axis=1)- np.nanmin(coordinates[:,2::3], axis=1))  
    
    return width
    
def calculate_boundary_quality(boundary_coordinates):
    """
    Calculate the standard deviation of the z-values of the 
    boundary cells. 
    
    """
    n_times = boundary_coordinates.shape[0]   
    
    z_values = boundary_coordinates[:, 2::3]
    n_cells_boundary = np.sum(~np.isnan(z_values[0,:]))#we assume that the number of cells in the boundary is constant over time   
    z_values_boundary = np.reshape(z_values[~np.isnan(z_values)], (n_times, n_cells_boundary))
    #mean_z_values = np.mean(z_values_boundary, axis=1)
    std_z_values = np.std(z_values_boundary, axis=1, ddof=1) 

    return std_z_values;
    
        


def calculate_column_quality(patch_width):
    """
    Calculate the column quality at a given time point by adding the width in 
    x and y direction and the deviation from the max height in z-direction.
    
    """
    #return patch_width[0] + patch_width[1] + abs(n_cells_in_patch - patch_width[2])
    return patch_width[:,0] + patch_width[:,1] 

def calculate_relative_column_quality(patch_coordinates):
    """
    Calculate the column quality at a given time point by adding the width in 
    x and y direction and the deviation from the max height in z-direction all 
    relative to the maximum width.    
    """    
    
    
    patch_width = calculate_width(patch_coordinates)
    n_cells = calculate_number_cells_in_patch(patch_coordinates)
    cell_diameter = 1.0
    target = (n_cells-1)*cell_diameter
    relative_column_quality = np.zeros(patch_coordinates.shape[0])
    relative_column_quality[n_cells != 1] = (target[n_cells != 1] - patch_width[(n_cells != 1),2] + patch_width[n_cells != 1,0]+ patch_width[n_cells != 1,1] )/(3*target[n_cells != 1])
    
    return relative_column_quality
    
    
def calculate_number_cells_in_patch(patch_coordinates):
    """
    Calculate the number of cells in a given patch over time from its coordinates.
    """
    
    return np.count_nonzero(~np.isnan(patch_coordinates[:,0::3]), axis=1) #don't forget to count only in one dimension.
    
def calculate_relative_column_height(patch_coordinates):
    """
    Calculates the relative column height of a given patch defined as the maximum
    difference in z values compared to n_cells*cell_diameter.
    """
    column_height = abs(np.nanmax(patch_coordinates[:,2::3], axis=1)- np.nanmin(patch_coordinates[:,2::3], axis=1)) 
    n_cells = calculate_number_cells_in_patch(patch_coordinates)
    
    cell_diameter = 1.0
    target_height = (n_cells-1)*cell_diameter #we subtract one in order to also only measure between midpoints 
    #target_height[target_height==0.0] =  1.0 # catch if the target height is zero (patches with a single cell)   
    
    relative_column_height = np.ones(patch_coordinates.shape[0])
    relative_column_height[n_cells != 1] = column_height[n_cells != 1]/target_height[n_cells != 1]      
    
    return relative_column_height
    


if __name__ == "__main__":
    
    path = 'results_from_time_0/'
    times, coordinates = read_in_data(path)
    print(coordinates.shape)
    print(coordinates.dtype)
    
    width = calculate_width(coordinates)
    plt.plot(times, width)
    plt.show()
    
    boundary_layers = get_border_layers(path)

    for boundary in boundary_layers:
        
        assert boundary.shape == coordinates.shape
        boundary_coordinates = coordinates.copy()
        boundary_coordinates[~boundary] = np.nan
        assert np.count_nonzero(~np.isnan(boundary_coordinates[0,:])) == 60
        
        #bw = calculate_width(boundary_coordinates)
        bs = calculate_boundary_quality(boundary_coordinates)
        lineObj = plt.plot(times, bs)
        plt.legend(lineObj, ('x', 'y', 'z'))
        plt.show()
#
#    clp = get_clonal_patches(path)
#    
#    
#    for ancestor, patch in clp.items():
#        print(ancestor)
#        assert patch.shape == coordinates.shape
#        if ancestor != -1:
#            pc = coordinates.copy()
#            pc[~patch] = np.nan
#        
#            pw = calculate_width(pc)
#            lineObj = plt.plot(times, pw)
#            plt.legend(lineObj, ('x', 'y', 'z'))
#            plt.show()
    