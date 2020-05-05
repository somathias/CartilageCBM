# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:17:38 2017

@author: kubuntu1404
"""

import multiprocessing
import os
import subprocess
import sys

#sys.path.append('/home/kubuntu1804/Documents/sheet_metrics_python')
sys.path.append('analysis')

import evaluate_cartilage_sheet

#import numpy as np

#executable = '/home/kubuntu1804/Documents/chaste_build/projects/cartilage/apps/CartilageSheetSimulation'
#output_path = '/home/kubuntu1804/Documents/sf_simulation_results/'



#number_of_simulations = 3

def main():
    output_directory = 'dev-cartilage_sheet_class/'
    run_simulations(output_directory)
    
    run_postprocessing(output_directory)
    
    print('Done.')


# Create a list of commands and pass them to separate processes
def run_simulations_multiple_random_seeds(output_directory, flags='',  number_of_simulations=3, executable='/home/kubuntu1804/Documents/chaste_build/projects/cartilage/apps/CartilageSheetSimulation'):
    
    if not(os.path.isfile(executable)):
        raise Exception('Could not find executable: ' + executable)

    # Make a list of calls to a Chaste executable
    command_list = []

    base_command = 'nice -n 19 ' + executable

    for random_seed in range(number_of_simulations):
        command = base_command + ' --output-dir ' + output_directory + ' --S ' + str(random_seed) + flags 
        command_list.append(command)
        print(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(execute_command, command_list).get(86400)
    
# Create a list of commands and pass them to separate processes
def run_simulations_list(list_of_flags=[''], random_seed=0, executable='/home/kubuntu1804/Documents/chaste_build/projects/cartilage/apps/CartilageSheetSimulation'):
    
    if not(os.path.isfile(executable)):
        raise Exception('Could not find executable: ' + executable)

    # Make a list of calls to a Chaste executable
    command_list = []

    base_command = 'nice -n 19 ' + executable

    for flags in list_of_flags:

        command = base_command + ' --S ' + str(random_seed) + flags 
        command_list.append(command)
        print(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(execute_command, command_list).get(86400)

def run_postprocessing_multiple_random_seeds(number_of_simulations, output_directory, output_path='/home/kubuntu1804/Documents/sf_simulation_results/', start_time=0):
    
    # Make a list of output_directories
    directory_list = []

    for random_seed in range(number_of_simulations):

        directory = output_path + output_directory + str(random_seed) +'/results_from_time_'+str(start_time)+'/'
        directory_list.append(directory)
    print(directory_list)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting postprocessing with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(evaluate_cartilage_sheet.main, directory_list).get(86400)

def run_postprocessing_list(list_of_directories, random_seed=0, output_path='/home/kubuntu1804/Documents/sf_simulation_results/', start_time=0):
    
    # Make a list of output_directories
    directory_list = []

    for direc in list_of_directories:
        directory_list.append(output_path + direc + str(random_seed) +'/results_from_time_'+str(start_time)+'/')
    print(directory_list)

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
