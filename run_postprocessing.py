# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:17:38 2017

@author: kubuntu1404
"""

import multiprocessing
import os
import subprocess
import sys
import time

sys.path.append('/home/kubuntu1404/Documents/sheet_metrics')

import evaluate_cartilage_sheet

#import numpy as np

executable = '/home/kubuntu1404/Documents/scaling_cartilage_sheets/apps/src/CartilageSheetSimulation'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

number_of_simulations = 5

def main():
    
    #spring_stiffness_values = [1.0, 5.0, 12.5, 25, 50]
    activation_percentage_values = [0.05, 0.25, 0.5, 0.75, 1.0]
    spring_stiffness_values = [0.01, 0.1, 0.5]
    #activation_percentage_values = [0.05]
    
    output_directory = 'OnlyPerichondrialLayers/'
    random_seed = 0   
    run_postprocessing('testoutput/'+output_directory, spring_stiffness_values, activation_percentage_values, random_seed)    
    random_seed = 1
    run_postprocessing('testoutput/'+output_directory, spring_stiffness_values, activation_percentage_values, random_seed) 
    random_seed = 2   
    run_postprocessing('testoutput/'+output_directory, spring_stiffness_values, activation_percentage_values, random_seed)    
    random_seed = 3
    run_postprocessing('testoutput/'+output_directory, spring_stiffness_values, activation_percentage_values, random_seed) 
    print('Done. '+time.strftime("%Y%m%d-%H%M%S"))



def run_postprocessing(output_directory, spring_stiffness_values, activation_percentage_values, random_seed):
    
    # Make a list of output_directories
    directory_list = []

    for spring_stiffness in spring_stiffness_values:
        
        for activation_percentage in activation_percentage_values:

            directory = output_directory + '/' + str(spring_stiffness) + '/' + str(activation_percentage) + '/' + str(random_seed) +'/results_from_time_0/'
            directory_list.append(directory)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting postprocessing with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(evaluate_cartilage_sheet.main, directory_list).get(86400)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()