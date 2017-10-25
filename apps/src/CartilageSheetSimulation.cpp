/*
 * CartilageSheetSimulation.cpp
 *
 *  Created on: Oct 25, 2017
 *      Author: Sonja Mathias
 */

#include <cxxtest/TestSuite.h>

// Includes from trunk
#include "CellId.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//#include "UniformlyDistributedCellCycleModel.hpp"

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

/*
 * Prototype functions
 */
void SetupSingletons(unsigned randomSeed);
void DestroySingletons();
void SetupAndRunSimulation(unsigned randomSeed);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a Chaste executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("S", boost::program_options::value<unsigned>()->default_value(0),"The random seed");

    // Define parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // Print help message if wanted
    if (variables_map.count("help"))
    {
        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // Get ID and name from command line
    unsigned random_seed = variables_map["S"].as<unsigned>();

    SetupSingletons(random_seed);
    SetupAndRunSimulation(random_seed);
    DestroySingletons();
}

void SetupSingletons(unsigned randomSeed)
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
    RandomNumberGenerator::Instance()->Reseed(randomSeed);
    CellPropertyRegistry::Instance()->Clear();
    CellId::ResetMaxCellId();
}

void DestroySingletons()
{
    // This is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
    CellPropertyRegistry::Instance()->Clear();
}

void SetupAndRunSimulation(unsigned randomSeed)
{
    // Simulation goes here

    //std::string output_directory = "paths/to/output" + boost::lexical_cast<std::string>(randomSeed);

    std::cout << "Testing2! " << randomSeed << std::endl;
}


