# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 13:06:55 2017

@author: Sonja Mathias
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import sheet_metrics as sm

plt.style.use('ggplot')


def main(path):
    ### General #################################################################
    
    #path = '/home/kubuntu1404/Documents/testoutput_171010_exp-essence_meeting/3dNodeBasedCartilageSheet/Test3dAdhesionVsSynchro/HighAdhesionSynchro/cutoff1-5/results_from_time_0/'
    #path = 'results_from_time_0/'
    
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'
    os.makedirs(results_path, exist_ok=True)
    
#    # calculate the sheet width and plot it as a function of time
#    width = sm.calculate_width(coordinates)
#    # normalize to initial width
#    initial_width = width[0,:] 
#    width_deviation = width - initial_width*np.ones(width.shape)
#    
#    # save to file
#    np.savetxt(results_path + 'data_complete_sheet.txt', width_deviation, fmt='%f')    
#    
#    plt.figure()
#    lineObjects = plt.plot(times, width_deviation)
#    plt.legend(lineObjects, ('x', 'y', 'z'))
#    plt.title('Deviation from initial width of complete cell sheet')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.03])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Width in cell diameters')
#    #plt.show()
#    plt.savefig(results_path + 'complete_sheet.png')
#    plt.savefig(results_path + 'complete_sheet.pdf')
#    plt.close()
    
    
#    ### Boundary ################################################################
#    boundary_layers = sm.get_border_layers(path)
#    
#    for index, boundary in enumerate(boundary_layers):
#        
#        boundary_layers = ['lower', 'upper']
#        
#        assert boundary.shape == coordinates.shape
#        boundary_coordinates = coordinates.copy()
#        boundary_coordinates[~boundary] = np.nan
#        #assert np.count_nonzero(~np.isnan(boundary_coordinates[0,:])) == 60
#        
#        bw = sm.calculate_width(boundary_coordinates)
#        # normalize to initial width
#        initial_width = bw[0,:] 
#        bw_deviation = bw - initial_width*np.ones(bw.shape)
#        # save to file
#        np.savetxt(results_path + 'data_'+ boundary_layers[index] + '_boundary.txt', bw_deviation, fmt='%f') 
#        
#        plt.figure()
#        lineObj = plt.plot(times, bw_deviation)
#        plt.legend(lineObj, ('x', 'y', 'z'))
#        plt.title('Deviation from initial width of '+ boundary_layers[index]+' boundary layer')
#        ax = plt.gca()
#        ttl = ax.title
#        ttl.set_position([.5, 1.03])
#        plt.xlabel('Time in hours')
#        plt.ylabel('Width in cell diameters')
#        #plt.show()
#        plt.savefig(results_path + boundary_layers[index] +'boundary.png')
#        plt.savefig(results_path + boundary_layers[index] +'boundary.pdf')
#        plt.close()
#        
#        bs = sm.calculate_boundary_quality(boundary_coordinates)
#        np.savetxt(results_path + 'data_'+ boundary_layers[index] + '_boundary_std.txt', bs, fmt='%f')
#        
#        plt.figure()
#        lineObj = plt.plot(times, bs)
#        plt.title('Standard deviation of the z values of the '+ boundary_layers[index]+' boundary layer')
#        ax = plt.gca()
#        ttl = ax.title
#        ttl.set_position([.5, 1.03])
#        plt.xlabel('Time in hours')
#        plt.ylabel('Width in cell diameters')
#        #plt.show()
#        plt.savefig(results_path + boundary_layers[index] +'boundary_std.png')
#        plt.savefig(results_path + boundary_layers[index] +'boundary_std.pdf')
#        plt.close()
#    
#    plt.close('all')
    
    ## Clonal Envelopes ########################################################
    clp = sm.get_clonal_patches(path)
    
    #print(clp[11.0])
    
    n_patches = len(clp)

    #store column quality metric
    cq = np.zeros((n_patches, times.shape[0]))
    #store number of cells
    nc = np.zeros((n_patches, times.shape[0]))
    #store relative column height
    ch = np.zeros((n_patches, times.shape[0]))
    # store relative column quality metric
    rcq = np.zeros((n_patches, times.shape[0]))
    # store patch projection area
    ppa = np.zeros((n_patches, times.shape[0]))
    ind = 0
    for ancestor, patch in clp.items():
        #print(ancestor)
        assert patch.shape == coordinates.shape
        if ancestor != -1:
            pc = coordinates.copy()
            pc[~patch] = np.nan
            pw = sm.calculate_width(pc)
            
            #calculate number of cells
            n_cells = sm.calculate_number_cells_in_patch(pc)
            nc[ind,:] = n_cells
            #calculate relative height
            #ch[ind,:] = sm.calculate_relative_column_height(pc)
            #calculate column quality
            #sc = sm.calculate_column_quality(pw)
            #cq[ind,:] = sc 
            #calculate relative column quality
            #rcq[ind,:] = sm.calculate_relative_column_quality(pc)   
            #calculate patch projection area
            ppa[ind, :] = pw[:,0] * pw[:,1]
            
            ind = ind+1
        else:
            cq = cq[:-1, :]
            nc = nc[:-1, :]
            ch = ch[:-1, :]
            rcq = rcq[:-1,:]
            ppa = ppa[:-1, :]
    print(ppa[0,:])
#            
#    #--calculate, plot and save mean column quality ----------------------------
#    cq_mean = np.mean(cq, axis=0)
#    cq_std = np.std(cq, axis=0)
#    # save to file
#    np.savetxt(results_path + 'data_column_quality_av.txt', cq_mean, fmt='%f')
#    np.savetxt(results_path + 'data_column_quality_std.txt', cq_std, fmt='%f')      
#    
#    plt.figure()
#    lineObj = plt.plot(times, cq_mean)
#    plt.fill_between(times, cq_mean-cq_std, cq_mean+cq_std, alpha=0.5,  edgecolor='#CC4F1B', facecolor='#FF9848')
#    plt.title('Average column quality over time (+- standard deviation)')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.03])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Column quality score in cell diameters')
#    plt.savefig(results_path + 'metric_average.png')
#    plt.savefig(results_path + 'metric_average.pdf')
#    plt.close('all')
#
#    
    #--calculate, plot and save average number of cells per clonal patch--------------
    nc_mean = np.nanmean(nc, axis=0)
    nc_std = np.nanstd(nc, axis=0)
    # save to file
    np.savetxt(results_path + 'data_ncells_av.txt', nc_mean, fmt='%f')
    np.savetxt(results_path + 'data_ncells_std.txt', nc_std, fmt='%f')      
#    
#    plt.figure()
#    lineObj = plt.plot(times, nc_mean)
#    plt.fill_between(times, nc_mean-nc_std, nc_mean+nc_std, alpha=0.5,   edgecolor='#CC4F1B', facecolor='#FF9848')
#    plt.title('Average number of cells per clonal patch over time \n (+- standard deviation)')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.01])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Number of cells')
#    plt.savefig(results_path + 'ncells_average.png')
#    plt.savefig(results_path + 'ncells_average.pdf')
#    plt.close('all')
#    
#    #--calculate, plot and save average relative column height --------------------
#    ch_mean = np.mean(ch, axis=0)
#    ch_std = np.std(ch, axis=0)
#    # save to file
#    np.savetxt(results_path + 'data_relative_column_height_av.txt', ch_mean, fmt='%f')
#    np.savetxt(results_path + 'data_relative_column_height_std.txt', ch_std, fmt='%f')      
#    
#    plt.figure()
#    lineObj = plt.plot(times, ch_mean)
#    plt.fill_between(times, ch_mean-ch_std, ch_mean+ch_std, alpha=0.5,   edgecolor='#CC4F1B', facecolor='#FF9848')
#    plt.title('Average relative height of clonal patch over time \n (+- standard deviation)')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.01])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Relative height')
#    plt.savefig(results_path + 'relative_column_height.png')
#    plt.savefig(results_path + 'relative_column_height.pdf')
#    plt.close('all')
#    
#    
#    #--calculate, plot and save average relative column quality ----------------
#    rcq_mean = np.mean(rcq, axis=0)
#    rcq_std = np.std(rcq, axis=0)
#    # save to file
#    np.savetxt(results_path + 'data_relative_column_quality_av.txt', rcq_mean, fmt='%f')
#    np.savetxt(results_path + 'data_relative_column_quality_std.txt', rcq_std, fmt='%f')      
#    
#    plt.figure()
#    lineObj = plt.plot(times, rcq_mean)
#    plt.fill_between(times, rcq_mean-rcq_std, rcq_mean+rcq_std, alpha=0.5,   edgecolor='#CC4F1B', facecolor='#FF9848')
#    plt.title('Average relative column quality over time \n (+- standard deviation)')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.01])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Relative column quality')
#    plt.savefig(results_path + 'relative_column_quality.png')
#    plt.savefig(results_path + 'relative_column_quality.pdf')
#    plt.close('all')
#    
    #--calculate and save mean patch projection area----------------------------
    ppa_mean = np.nanmean(ppa, axis=0)
    ppa_std = np.nanstd(ppa, axis=0)
    # save to file
    np.savetxt(results_path + 'data_patch_projection_area_av.txt', ppa_mean, fmt='%f')
    np.savetxt(results_path + 'data_patch_projection_area_std.txt', ppa_std, fmt='%f')  
#    
    #-- calculate, plot and save patch size distribution --------------------------------------
    # find the values at t=35h and t=40h
    ind_t35 = np.where(times==35.0)[0][0] #what if we did not sample at this time???
    ind_t40 = np.where(times==40.0)[0][0]
    ind_t45 = np.where(times==45.0)[0][0]

    max_cells_per_all_columns = int(nc[:, ind_t45].max())
    bins=np.arange(1,max_cells_per_all_columns+2) #we add two because the last bin has to include max and max+1, else it would contain max-1 and max
    
    plt.figure()
    plt.hist([nc[:, ind_t35], nc[:, ind_t40], nc[:, ind_t45]], bins=bins, align='left', color= ['#fcbba1', '#fb6a4a', '#a50f15'])
    ax = plt.gca()
    ax.set_xticks(bins)
    ax.set_xticklabels(['1', '2', '3','4','5','6','7','8','9','10','11'])
    plt.legend(('t=35h', 't=40h', 't=45h'))
    plt.xlabel('Size of clonal patch in number of cells')
    plt.ylabel('Number of clonal patches')
    plt.savefig(results_path + 'patch_size_distribution.png')
    plt.savefig(results_path + 'patch_size_distribution.pdf')

    
#    patch_size_dist_t35, bin_edges_t35 = np.histogram(nc[:, ind_t35], bins=bins)
#    patch_size_dist_t40, bin_edges_t40 = np.histogram(nc[:, ind_t40], bins=bins) 
#    patch_size_dist_t45, bin_edges_t45 = np.histogram(nc[:, ind_t45], bins=bins)
#
#    # make a boolean mask for each of the cell counts
#    len_counts=len(bins[:-1])
#    rch_t35_mean_per_cell_count = np.zeros(len_counts)
#    rch_t40_mean_per_cell_count = np.zeros(len_counts)
#    rch_t45_mean_per_cell_count = np.zeros(len_counts)
#    rcq_t35_mean_per_cell_count = np.zeros(len_counts)
#    rcq_t40_mean_per_cell_count = np.zeros(len_counts)
#    rcq_t45_mean_per_cell_count = np.zeros(len_counts)
#
#    for count in bins[:-1]:
#        #t=35h
#        mask = nc[:, ind_t35] == count  
#        rch_t35 = ch[:, ind_t35]
#        rch_t35_mean_per_cell_count[count-1] = np.mean(rch_t35[mask])
#        rcq_t35 = rcq[:, ind_t35]
#        rcq_t35_mean_per_cell_count[count-1] = np.mean(rcq_t35[mask])
#        
#        #t=40h
#        mask = nc[:, ind_t40] == count  
#        rch_t40 = ch[:, ind_t40]
#        rch_t40_mean_per_cell_count[count-1] = np.mean(rch_t40[mask])
#        rcq_t40 = rcq[:, ind_t40]
#        rcq_t40_mean_per_cell_count[count-1] = np.mean(rcq_t40[mask])
#        
#        #t=45h
#        mask = nc[:, ind_t45] == count 
#        rch_t45 = ch[:, ind_t45]
#        rch_t45_mean_per_cell_count[count-1] = np.mean(rch_t45[mask])
#        rcq_t45 = rcq[:, ind_t45]
#        rcq_t45_mean_per_cell_count[count-1] = np.mean(rcq_t45[mask])
#   
#    data_t35 = np.array([bins[:-1], patch_size_dist_t35, rch_t35_mean_per_cell_count, rcq_t35_mean_per_cell_count]).T
#    data_t40 = np.array([bins[:-1], patch_size_dist_t40, rch_t40_mean_per_cell_count, rcq_t40_mean_per_cell_count]).T
#    data_t45 = np.array([bins[:-1], patch_size_dist_t45, rch_t45_mean_per_cell_count, rcq_t45_mean_per_cell_count]).T
#    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t35.txt', data_t35, fmt='%f')
#    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t40.txt', data_t40, fmt='%f')
#    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t45.txt', data_t45, fmt='%f')
#
#    plt.figure()
##    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, bins[:-1], rcq_t40_mean_per_cell_count, bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= ['#fcbba1', '#fb6a4a', '#a50f15'])
#    plt.plot(bins[:-1], rch_t35_mean_per_cell_count, linestyle='None', marker='o',  color= '#fcbba1')
#    plt.plot(bins[:-1], rch_t40_mean_per_cell_count, linestyle='None', marker='o',  color= '#fb6a4a')
#    plt.plot(bins[:-1], rch_t45_mean_per_cell_count, linestyle='None', marker='o',  color= '#a50f15')
#    plt.legend(('t=35h', 't=40h', 't=45h'))
#    plt.xlabel('Size of clonal patch in number of cells')
#    plt.ylabel('Average relative column height')
#    plt.savefig(results_path + 'rel_column_height_as_function_of_patch_size.png')
#    plt.savefig(results_path + 'rel_column_height_as_function_of_patch_size.pdf')
#    
#    
#    plt.figure()
##    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, bins[:-1], rcq_t40_mean_per_cell_count, bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= ['#fcbba1', '#fb6a4a', '#a50f15'])
#    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, linestyle='None', marker='o',  color= '#fcbba1')
#    plt.plot(bins[:-1], rcq_t40_mean_per_cell_count, linestyle='None', marker='o',  color= '#fb6a4a')
#    plt.plot(bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= '#a50f15')
#    plt.legend(('t=35h', 't=40h', 't=45h'))
#    plt.xlabel('Size of clonal patch in number of cells')
#    plt.ylabel('Average relative column quality')
#    plt.savefig(results_path + 'rel_column_quality_as_function_of_patch_size.png')
#    plt.savefig(results_path + 'rel_column_quality_as_function_of_patch_size.pdf')
#    #plt.show()
    plt.close('all')
    

def update_column_quality_data(path):
       
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'
    os.makedirs(results_path, exist_ok=True)
    
    ## Clonal Envelopes ########################################################
    clp = sm.get_clonal_patches(path)
    
    n_patches = len(clp)

    #store column quality metric
    cq = np.zeros((n_patches, times.shape[0]))
    ind = 0
    for ancestor, patch in clp.items():
        #print(ancestor)
        assert patch.shape == coordinates.shape
        if ancestor != -1:
            pc = coordinates.copy()
            pc[~patch] = np.nan
            pw = sm.calculate_width(pc)
            
            #calculate column quality
            sc = sm.calculate_column_quality(pw)
            cq[ind,:] = sc                    
            ind = ind+1
        else:
            cq = cq[:-1, :]

    # save to file
    np.savetxt(results_path + 'data_column_quality_all.txt', cq, fmt='%f')
     
    print('Test')
    plt.figure()
    plt.plot(times, cq.T)
    #plt.show()
#    plt.fill_between(times, cq_mean-cq_std, cq_mean+cq_std, alpha=0.5,  edgecolor='#CC4F1B', facecolor='#FF9848')
#    plt.title('Average column quality over time (+- standard deviation)')
#    ax = plt.gca()
#    ttl = ax.title
#    ttl.set_position([.5, 1.03])
#    plt.xlabel('Time in hours')
#    plt.ylabel('Column quality score in cell diameters')
#    plt.savefig(results_path + 'metric_average.png')
#    plt.savefig(results_path + 'metric_average.pdf')
    plt.close('all')
    
def update_boundary_metric(path):
    
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'
    os.makedirs(results_path, exist_ok=True)
    
    ### Boundary ################################################################
    boundary_layers = sm.get_border_layers(path)
    
    for index, boundary in enumerate(boundary_layers):
        
        boundary_layers = ['lower', 'upper']
        
        assert boundary.shape == coordinates.shape
        boundary_coordinates = coordinates.copy()
        boundary_coordinates[~boundary] = np.nan
        #assert np.count_nonzero(~np.isnan(boundary_coordinates[0,:])) == 60
        
        bs = sm.calculate_boundary_quality(boundary_coordinates)
        np.savetxt(results_path + 'data_'+ boundary_layers[index] + '_boundary_std.txt', bs, fmt='%f')
        
        plt.figure()
        plt.plot(times, bs)
        plt.title('Standard deviation of the z values of the '+ boundary_layers[index]+' boundary layer')
        ax = plt.gca()
        ttl = ax.title
        ttl.set_position([.5, 1.03])
        plt.xlabel('Time in hours')
        plt.ylabel('Width in cell diameters')
        #plt.show()
        plt.savefig(results_path + boundary_layers[index] +'boundary_std.png')
        plt.savefig(results_path + boundary_layers[index] +'boundary_std.pdf')
        plt.close()
    
    plt.close('all')
    
def update_relative_column_height(path):
    ## Clonal Envelopes ########################################################
    
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'

    clp = sm.get_clonal_patches(path)
    
    n_patches = len(clp)

    #store relative column height
    ch = np.zeros((n_patches, times.shape[0]))
    ind = 0
    for ancestor, patch in clp.items():
        #print(ancestor)
        assert patch.shape == coordinates.shape
        if ancestor != -1:
            pc = coordinates.copy()
            pc[~patch] = np.nan

            #calculate relative height
            ch[ind,:] = sm.calculate_relative_column_height(pc)          
            
            ind = ind+1
        else:
            ch = ch[:-1, :]
    
    #calculate, plot and save average relative column height
    ch_mean = np.mean(ch, axis=0)
    ch_std = np.std(ch, axis=0)
    # save to file
    np.savetxt(results_path + 'data_relative_column_height_av.txt', ch_mean, fmt='%f')
    np.savetxt(results_path + 'data_relative_column_height_std.txt', ch_std, fmt='%f')      
    
    plt.figure()
    plt.plot(times, ch_mean)
    plt.fill_between(times, ch_mean-ch_std, ch_mean+ch_std, alpha=0.5,   edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.title('Average relative height of clonal patch over time \n (+- standard deviation)')
    ax = plt.gca()
    ttl = ax.title
    ttl.set_position([.5, 1.01])
    plt.xlabel('Time in hours')
    plt.ylabel('Relative height')
    plt.savefig(results_path + 'relative_column_height.png')
    plt.savefig(results_path + 'relative_column_height.pdf')
    plt.close('all')    
    
def update_relative_column_quality(path):
    ## Clonal Envelopes ########################################################
    
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'

    clp = sm.get_clonal_patches(path)
    
    n_patches = len(clp)

    #store relative column quality
    rcq = np.zeros((n_patches, times.shape[0]))
    ind = 0
    for ancestor, patch in clp.items():
        #print(ancestor)
        assert patch.shape == coordinates.shape
        if ancestor != -1:
            pc = coordinates.copy()
            pc[~patch] = np.nan

            #calculate relative height
            rcq[ind,:] = sm.calculate_relative_column_quality(pc)          
            
            ind = ind+1
        else:
            rcq = rcq[:-1, :]
    
    #calculate, plot and save average relative column quality
    rcq_mean = np.mean(rcq, axis=0)
    rcq_std = np.std(rcq, axis=0)
    # save to file
    np.savetxt(results_path + 'data_relative_column_quality_av.txt', rcq_mean, fmt='%f')
    np.savetxt(results_path + 'data_relative_column_quality_std.txt', rcq_std, fmt='%f')      
    
    plt.figure()
    plt.plot(times, rcq_mean)
    plt.fill_between(times, rcq_mean-rcq_std, rcq_mean+rcq_std, alpha=0.5,   edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.title('Average relative column quality over time \n (+- standard deviation)')
    ax = plt.gca()
    ttl = ax.title
    ttl.set_position([.5, 1.01])
    plt.xlabel('Time in hours')
    plt.ylabel('Relative column quality')
    plt.savefig(results_path + 'relative_column_quality.png')
    plt.savefig(results_path + 'relative_column_quality.pdf')
    plt.close('all')
    
def update_patch_size_distribution(path):
    ## Clonal Envelopes ########################################################
    
    # read in the data 
    times, coordinates = sm.read_in_data(path);
    
    
    results_path = path + 'metrics_results/'

    clp = sm.get_clonal_patches(path)
    
    n_patches = len(clp)
    
    #store number of cells
    nc = np.zeros((n_patches, times.shape[0]))
    #store relative column height
    ch = np.zeros((n_patches, times.shape[0]))
    #store relative column quality
    rcq = np.zeros((n_patches, times.shape[0]))
    ind = 0
    for ancestor, patch in clp.items():
        #print(ancestor)
        assert patch.shape == coordinates.shape
        if ancestor != -1:
            pc = coordinates.copy()
            pc[~patch] = np.nan
                         
            #calculate number of cells
            n_cells = sm.calculate_number_cells_in_patch(pc)
            nc[ind,:] = n_cells
            #calculate relative height
            ch[ind,:] = sm.calculate_relative_column_height(pc)
            #calculate relative height
            rcq[ind,:] = sm.calculate_relative_column_quality(pc)          
            
            ind = ind+1
        else:
            nc = nc[:-1, :]
            ch = ch[:-1, :]
            rcq = rcq[:-1, :]
            
    # calculate patch size distribution
    # find the values at t=35h and t=40h
    ind_t35 = np.where(times==35.0)[0][0] #what if we did not sample at this time???
    ind_t40 = np.where(times==40.0)[0][0]
    ind_t45 = np.where(times==45.0)[0][0]

    max_cells_per_all_columns = int(nc[:, ind_t45].max())
    bins=np.arange(1,max_cells_per_all_columns+2) #we add two because the last bin has to include max and max+1, else it would contain max-1 and max
    
    plt.figure()
    plt.hist([nc[:, ind_t35], nc[:, ind_t40], nc[:, ind_t45]], bins=bins, align='left', color= ['#fcbba1', '#fb6a4a', '#a50f15'])
    ax = plt.gca()
    ax.set_xticks(bins)
    ax.set_xticklabels(['1', '2', '3','4','5','6','7','8','9','10','11'])
    plt.legend(('t=35h', 't=40h', 't=45h'))
    plt.xlabel('Size of clonal patch in number of cells')
    plt.ylabel('Number of clonal patches')
    plt.savefig(results_path + 'patch_size_distribution.png')
    plt.savefig(results_path + 'patch_size_distribution.pdf')
    
    patch_size_dist_t35, bin_edges_t35 = np.histogram(nc[:, ind_t35], bins=bins)
    patch_size_dist_t40, bin_edges_t40 = np.histogram(nc[:, ind_t40], bins=bins) 
    patch_size_dist_t45, bin_edges_t45 = np.histogram(nc[:, ind_t45], bins=bins)

    # make a boolean mask for each of the cell counts
    len_counts=len(bins[:-1])
    rch_t35_mean_per_cell_count = np.zeros(len_counts)
    rch_t40_mean_per_cell_count = np.zeros(len_counts)
    rch_t45_mean_per_cell_count = np.zeros(len_counts)
    rcq_t35_mean_per_cell_count = np.zeros(len_counts)
    rcq_t40_mean_per_cell_count = np.zeros(len_counts)
    rcq_t45_mean_per_cell_count = np.zeros(len_counts)

    for count in bins[:-1]:
        #t=35h
        mask = nc[:, ind_t35] == count  
        rch_t35 = ch[:, ind_t35]
        rch_t35_mean_per_cell_count[count-1] = np.mean(rch_t35[mask])
        rcq_t35 = rcq[:, ind_t35]
        rcq_t35_mean_per_cell_count[count-1] = np.mean(rcq_t35[mask])
        
        #t=40h
        mask = nc[:, ind_t40] == count  
        rch_t40 = ch[:, ind_t40]
        rch_t40_mean_per_cell_count[count-1] = np.mean(rch_t40[mask])
        rcq_t40 = rcq[:, ind_t40]
        rcq_t40_mean_per_cell_count[count-1] = np.mean(rcq_t40[mask])
        
        #t=45h
        mask = nc[:, ind_t45] == count 
        rch_t45 = ch[:, ind_t45]
        rch_t45_mean_per_cell_count[count-1] = np.mean(rch_t45[mask])
        rcq_t45 = rcq[:, ind_t45]
        rcq_t45_mean_per_cell_count[count-1] = np.mean(rcq_t45[mask])
               
               
    data_t35 = np.array([bins[:-1], patch_size_dist_t35, rch_t35_mean_per_cell_count, rcq_t35_mean_per_cell_count]).T
    data_t40 = np.array([bins[:-1], patch_size_dist_t40, rch_t40_mean_per_cell_count, rcq_t40_mean_per_cell_count]).T
    data_t45 = np.array([bins[:-1], patch_size_dist_t45, rch_t45_mean_per_cell_count, rcq_t45_mean_per_cell_count]).T
    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t35.txt', data_t35, fmt='%f')
    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t40.txt', data_t40, fmt='%f')
    np.savetxt(results_path + 'data_column_quality_as_function_of_patch_size_t45.txt', data_t45, fmt='%f')
    

    plt.figure()
#    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, bins[:-1], rcq_t40_mean_per_cell_count, bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= ['#fcbba1', '#fb6a4a', '#a50f15'])
    plt.plot(bins[:-1], rch_t35_mean_per_cell_count, linestyle='None', marker='o',  color= '#fcbba1')
    plt.plot(bins[:-1], rch_t40_mean_per_cell_count, linestyle='None', marker='o',  color= '#fb6a4a')
    plt.plot(bins[:-1], rch_t45_mean_per_cell_count, linestyle='None', marker='o',  color= '#a50f15')
    plt.legend(('t=35h', 't=40h', 't=45h'))
    plt.xlabel('Size of clonal patch in number of cells')
    plt.ylabel('Average relative column height')
    plt.savefig(results_path + 'rel_column_height_as_function_of_patch_size.png')
    plt.savefig(results_path + 'rel_column_height_as_function_of_patch_size.pdf')
    
    plt.figure()
#    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, bins[:-1], rcq_t40_mean_per_cell_count, bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= ['#fcbba1', '#fb6a4a', '#a50f15'])
    plt.plot(bins[:-1], rcq_t35_mean_per_cell_count, linestyle='None', marker='o',  color= '#fcbba1')
    plt.plot(bins[:-1], rcq_t40_mean_per_cell_count, linestyle='None', marker='o',  color= '#fb6a4a')
    plt.plot(bins[:-1], rcq_t45_mean_per_cell_count, linestyle='None', marker='o',  color= '#a50f15')
    plt.legend(('t=35h', 't=40h', 't=45h'))
    plt.xlabel('Size of clonal patch in number of cells')
    plt.ylabel('Average relative column quality')
    plt.savefig(results_path + 'rel_column_quality_as_function_of_patch_size.png')
    plt.savefig(results_path + 'rel_column_quality_as_function_of_patch_size.pdf')


    
if __name__ == "__main__":
    #main('test_data/results_from_time_0/')
    main('/home/kubuntu1804/Documents/sf_simulation_results/dev/0/results_from_time_0/')