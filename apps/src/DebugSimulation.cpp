/*
 * DebugSimulation.cpp
 *
 *  Created on: Oct 25, 2017
 *      Author: Sonja Mathias
 */

#include <cxxtest/TestSuite.h>
#include <string>
#include <iostream>

// Includes from trunk
#include "CellId.hpp"
//#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
//#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
//#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//#include "UniformlyDistributedCellCycleModel.hpp"

#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "NodeBasedCartilageSheet.hpp"
#include "PatchSizeTrackingModifier.hpp"


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

int main(int argc, char *argv[]) {
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::InitializePetsc(&argc, &argv);

    unsigned random_seed = 67;
    SetupSingletons(random_seed);

	// Construct a new cartilage sheet
    NodeBasedCartilageSheet* p_cartilage_sheet =
            new NodeBasedCartilageSheet();

    // set the sheet dimensions
    p_cartilage_sheet->SetCartilageSheetDimensions(5, 5, 2);

    // generate the nodes
    p_cartilage_sheet->GenerateNodesOnStackedHexagonalGrid(1.0);

    // setup the cell population
    // if (!p_cartilage_sheet->isCellPopulationSetup()) {
    // 	p_cartilage_sheet->Setup();
    // }
    p_cartilage_sheet->Setup();

    // Delete middle column
    p_cartilage_sheet->InitialiseMissingColumnExperiment();
    std::cout << "Passed initialization." << std::endl;

    // get the cell population
    boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
    p_cartilage_sheet->GetCellPopulation();

    //pass cell population to the simulator
    OffLatticeSimulation<3> simulator(*cell_population);
    simulator.SetOutputDirectory("TestMissingColumnExperiment");
    simulator.SetSamplingTimestepMultiple(12);
    simulator.SetEndTime(50.0);

    MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
    p_force->SetCutOffLength(1.5);
    p_force->SetMeinekeSpringStiffness(1.0);
    simulator.AddForce(p_force);

    // Add the PatchSizeTracker to ensure that patches have maximum 6 cells
    // Currently it is not possible to run simulations with the CellTissueTypeBasedCellCycleModel in them without this.
    MAKE_PTR(PatchSizeTrackingModifier<3>, p_modifier);
    simulator.AddSimulationModifier(p_modifier);

    simulator.Solve();

	DestroySingletons();
}

void SetupSingletons(unsigned randomSeed) {
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);

	// Reseed with 0 for same random numbers each time, or time(NULL) or simulation_id to change each realisation
	RandomNumberGenerator::Instance()->Reseed(randomSeed);
	CellPropertyRegistry::Instance()->Clear();
	CellId::ResetMaxCellId();
}

void DestroySingletons() {
	// This is from the tearDown method of the test suite
	SimulationTime::Destroy();
	RandomNumberGenerator::Destroy();
	CellPropertyRegistry::Instance()->Clear();
}




