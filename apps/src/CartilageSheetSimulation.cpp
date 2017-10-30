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
void SetupAndRunCartilageSheetSimulation(unsigned randomSeed, unsigned,
		unsigned, unsigned, double, double, double, std::string);

int main(int argc, char *argv[]) {
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::InitializePetsc(&argc, &argv);


	// Define command line options
	boost::program_options::options_description general_options(
			"This is a Chaste executable.\n");
	general_options.add_options()("help", "Produce help message")("S",
			boost::program_options::value<unsigned>()->default_value(0),
			"The random seed")("sw",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in x direction")("sd",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in y direction")("sh",
			boost::program_options::value<unsigned>()->default_value(2),
			"The number of cells in z direction")("mu",
			boost::program_options::value<double>()->default_value(15.0),
			"The spring stiffness")("A",
			boost::program_options::value<double>()->default_value(0.5),
			"The percentage of activated stem cells")("T",
			boost::program_options::value<double>()->default_value(10.0),
			"The simulation end time")("output-dir",
			boost::program_options::value<std::string>()->default_value(
					"3dNodeBasedCartilageSheet/"), "The output directory");

	// Define parse command line into variables_map
	boost::program_options::variables_map variables_map;
	boost::program_options::store(
			parse_command_line(argc, argv, general_options), variables_map);

	// Print help message if wanted
	if (variables_map.count("help")) {
		//std::cout << setprecision(3) << general_options << "\n";
		std::cout << general_options << "\n";
		return 1;
	}

	// Get ID and name from command line
	unsigned random_seed = variables_map["S"].as<unsigned>();
	unsigned n_cells_wide = variables_map["sw"].as<unsigned>();
	unsigned n_cells_deep = variables_map["sd"].as<unsigned>();
	unsigned n_cells_high = variables_map["sh"].as<unsigned>();
	double activation_percentage = variables_map["A"].as<double>();
	double spring_stiffness = variables_map["mu"].as<double>();
	double simulation_end_time = variables_map["T"].as<double>();
	std::string output_directory = variables_map["output-dir"].as<std::string>();

	SetupSingletons(random_seed);
	SetupAndRunCartilageSheetSimulation(random_seed, n_cells_wide, n_cells_deep,
			n_cells_high, activation_percentage, spring_stiffness,
			simulation_end_time, output_directory);
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

void SetupAndRunCartilageSheetSimulation(unsigned random_seed,
		unsigned n_cells_wide, unsigned n_cells_deep, unsigned n_cells_high,
		double activation_percentage, double spring_stiffness,
		double simulation_endtime, std::string output_directory) {
	// Simulation goes here

	bool random_birth_times = true;
	output_directory.append(boost::lexical_cast<std::string>(random_seed));

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
	unsigned n_activated_stem_cells = floor(
			n_cells_wide * n_cells_deep * activation_percentage);
	p_cartilage_sheet->InitialiseRandomStemCellConfiguration(
			n_activated_stem_cells);

	// get the cell population
	boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
			p_cartilage_sheet->GetCellPopulation();

	OffLatticeSimulationDirectedDivision<3> simulator(*cell_population);
	//OffLatticeSimulation<3> simulator(cell_population);
	simulator.SetOutputDirectory(output_directory);
	simulator.SetEndTime(simulation_endtime); //hours
	simulator.SetSamplingTimestepMultiple(12);

	MAKE_PTR(CellTissueTypeBasedGeneralisedLinearSpringForce<3>, p_force);
	double alpha = -2.0 * log(2.0 / spring_stiffness * 0.001);
	p_force->SetCutOffLength(1.5);
	p_force->SetMeinekeSpringStiffness(spring_stiffness);
	p_force->SetAlpha(alpha);
	simulator.AddForce(p_force);

	CellBasedEventHandler::Reset();
	simulator.Solve();

	CellBasedEventHandler::Headings();
	CellBasedEventHandler::Report();
}
