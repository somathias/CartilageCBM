/*
 * CartilageSheetSimulation.cpp
 *
 *  Created on: Oct 25, 2017
 *      Author: Sonja Mathias
 */

#include <cxxtest/TestSuite.h>
#include <string>

// Includes from trunk
#include "CellId.hpp"
//#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
//#include "DifferentiatedCellProliferativeType.hpp"
#include "ExecutableSupport.hpp"
//#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//#include "UniformlyDistributedCellCycleModel.hpp"

#include "CellBasedEventHandler.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "NodeBasedCartilageSheet.hpp"
#include "OffLatticeSimulationDirectedDivision.hpp"
#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"

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
void SetupAndRunCartilageSheetSimulation(unsigned randomSeed);

int main(int argc, char *argv[]) {
	// This sets up PETSc and prints out copyright information, etc.
	ExecutableSupport::StandardStartup(&argc, &argv);

	// Define command line options
	boost::program_options::options_description general_options(
			"This is a Chaste executable.\n");
	general_options.add_options()("help", "produce help message")("S",
			boost::program_options::value<unsigned>()->default_value(0),
			"The random seed");

	// Define parse command line into variables_map
	boost::program_options::variables_map variables_map;
	boost::program_options::store(
			parse_command_line(argc, argv, general_options), variables_map);

	// Print help message if wanted
	if (variables_map.count("help")) {
		std::cout << setprecision(3) << general_options << "\n";
		std::cout << general_options << "\n";
		return 1;
	}

	// Get ID and name from command line
	unsigned random_seed = variables_map["S"].as<unsigned>();

	SetupSingletons(random_seed);
	SetupAndRunCartilageSheetSimulation(random_seed);
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

void SetupAndRunCartilageSheetSimulation(unsigned random_seed) {
	// Simulation goes here

	//std::string output_directory = "paths/to/output" + boost::lexical_cast<std::string>(randomSeed);

	std::cout << "Testing3! " << random_seed << std::endl;

	unsigned n_cells_wide = 3;
	unsigned n_cells_deep = 3;
	unsigned n_cells_high = 2;
	bool random_birth_times = true;
	double spring_stiffness = 15;
	double activation_percentage =0.5;
	double alpha = 5.0 ;
	std::string output_directory = "3dNodeBasedCartilageSheet/" + boost::lexical_cast<std::string>(random_seed);;
	double simulation_endtime = 10.0;

	CellBasedEventHandler::Enable();

	NodeBasedCartilageSheet* p_cartilage_sheet = new NodeBasedCartilageSheet();

	// set the sheet dimensions
	p_cartilage_sheet->SetCartilageSheetDimensions(n_cells_wide, n_cells_deep,
			n_cells_high);
	p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);

	if (!random_birth_times) {
		p_cartilage_sheet->setSynchronizeCellCycles(true);
	}
	// generate the nodes
	p_cartilage_sheet->GenerateNodesOnHCPGrid();

	// setup the cell population
	if (!p_cartilage_sheet->isCellPopulationSetup()) {
		p_cartilage_sheet->Setup();
	}
	// setup the cell tissue types and cell division directions
	p_cartilage_sheet->InitialiseTissueLayersAndCellDivisionDirections();
	// setup the initial stem cell configuration
	unsigned n_activated_stem_cells = floor(n_cells_wide*n_cells_deep*activation_percentage);
	p_cartilage_sheet->InitialiseRandomStemCellConfiguration(n_activated_stem_cells);

	// get the cell population
	boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
			p_cartilage_sheet->GetCellPopulation();

	OffLatticeSimulationDirectedDivision<3> simulator(*cell_population);
	//OffLatticeSimulation<3> simulator(cell_population);
	simulator.SetOutputDirectory(output_directory);
	simulator.SetEndTime(simulation_endtime); //hours
	simulator.SetSamplingTimestepMultiple(12);

	MAKE_PTR(CellTissueTypeBasedGeneralisedLinearSpringForce<3>, p_force);
	p_force->SetCutOffLength(1.5);
	p_force->SetMeinekeSpringStiffness(spring_stiffness);
	p_force->SetAlpha(alpha);
	simulator.AddForce(p_force);

	CellBasedEventHandler::Reset();
	simulator.Solve();

	CellBasedEventHandler::Headings();
	CellBasedEventHandler::Report();
}

