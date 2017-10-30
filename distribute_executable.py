# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 09:17:38 2017

@author: kubuntu1404
"""

import multiprocessing
import os
import subprocess

import numpy as np

executable = '/home/kubuntu1404/Documents/scaling_cartilage_sheets/apps/src/CartilageSheetSimulation'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

number_of_simulations = 5

def main():
    output_directory = '3dNodeBasedCartilageSheet/'
    run_simulations(output_directory)
    run_postprocessing(output_directory)


# Create a list of commands and pass them to separate processes
def run_simulations(output_directory):

    # Make a list of calls to a Chaste executable
    command_list = []

    base_command = 'nice -n 19 ' + executable

    for random_seed in range(number_of_simulations):

        command = base_command + ' --output-dir ' + output_directory + ' --S ' + str(random_seed) 
        command_list.append(command)

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(execute_command, command_list).get(86400)

def run_postprocessing(output_directory):
    print("Test")


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()