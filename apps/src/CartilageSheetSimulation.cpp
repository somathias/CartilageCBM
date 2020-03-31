/*
 * CartilageSheetSimulation.cpp
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

#include "CellBasedEventHandler.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBoundaryCondition.hpp"

#include "NodeBasedCartilageSheet.hpp"
#include "OffLatticeSimulationDirectedDivision.hpp"
#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
#include "IndividualSpringStiffnessGeneralisedLinearSpringForce.hpp"
#include "CubicGeneralisedLinearSpringForce.hpp"
#include "RepulsionCubicForce.hpp"
#include "PWQGeneralisedLinearSpringForce.hpp"
#include "RepulsionForce.hpp"
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
void SetupAndRunCartilageSheetSimulation(unsigned randomSeed, bool, bool, bool, unsigned,
		unsigned, unsigned, unsigned, unsigned, unsigned, double, double, double, double, double, double, double, double,
		std::string, std::string);
void SetForceFunction(OffLatticeSimulation<3>&, std::string, double, double, double, double, double, double);

int main(int argc, char *argv[]) {
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::InitializePetsc(&argc, &argv);

	// Define command line options
	boost::program_options::options_description general_options(
			"This is a Chaste executable.\n");
	general_options.add_options()("help", "Produce help message")("sbt",
			"Synchronized birth times")("rdd",
			"Random division directions")("c",
			"Cartesian grid")("S",
			boost::program_options::value<unsigned>()->default_value(0),
			"The random seed")("sw",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in x direction")("sd",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in y direction")("sh",
			boost::program_options::value<unsigned>()->default_value(2),
			"The number of cells in z direction")("pl",
			boost::program_options::value<unsigned>()->default_value(1),
			"The number of lower perichondrial layers")("pu",
			boost::program_options::value<unsigned>()->default_value(0),
			"The number of upper perichondrial layers")("nb",
			boost::program_options::value<unsigned>()->default_value(0),
			"The number of rigid boundaries (0,1 or 2)")("u",
			boost::program_options::value<double>()->default_value(7.0),
			"The distance of the upper boundary to the lower one")("mu",
			boost::program_options::value<double>()->default_value(15.0),
			"The adhesion spring stiffness")("mu_R",
			boost::program_options::value<double>()->default_value(1.4),
			"The repulsion spring stiffness")("c",
			boost::program_options::value<double>()->default_value(1.0),
			"The homotypic chondrocyte multiplier")("b",
			boost::program_options::value<double>()->default_value(0.1),
			"The baseline directional adhesion multiplier")("A",
			boost::program_options::value<double>()->default_value(0.5),
			"The percentage of activated stem cells")("p",
			boost::program_options::value<double>()->default_value(0.0),
			"The maximum perturbation of the initial coordinates.")("T",
			boost::program_options::value<double>()->default_value(10.0),
			"The simulation end time")("F",
			boost::program_options::value<std::string>()->default_value(
					"CellTissueTypeBasedGLS"), "The force function used")("output-dir",
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

	bool random_birth_times = true;
	// set random birth times to false = synchronized case if wanted
	if (variables_map.count("sbt")) {
		random_birth_times = false;
	}
	bool random_division_directions = false;
	if (variables_map.count("rdd")) {
		random_division_directions = true;
	}
	bool cartesian_grid = false;
	if (variables_map.count("c")) {
		cartesian_grid = true;
	}

	// Get ID and name from command line
	unsigned random_seed = variables_map["S"].as<unsigned>();
	unsigned n_cells_wide = variables_map["sw"].as<unsigned>();
	unsigned n_cells_deep = variables_map["sd"].as<unsigned>();
	unsigned n_cells_high = variables_map["sh"].as<unsigned>();
	unsigned n_peri_lower = variables_map["pl"].as<unsigned>();
	unsigned n_peri_upper = variables_map["pu"].as<unsigned>();
	unsigned n_boundaries = variables_map["nb"].as<unsigned>();
	double upper_boundary = variables_map["u"].as<double>();
	double activation_percentage = variables_map["A"].as<double>();
	double maximum_perturbation = variables_map["p"].as<double>();
	double spring_stiffness = variables_map["mu"].as<double>();
	double spring_stiffness_repulsion = variables_map["mu_R"].as<double>();
	double homotypic_chondro_multiplier = variables_map["c"].as<double>();
	double baseline_adhesion_multiplier = variables_map["b"].as<double>();
	double simulation_end_time = variables_map["T"].as<double>();
	std::string force_function =
			variables_map["F"].as<std::string>();
	std::string output_directory =
			variables_map["output-dir"].as<std::string>();


	SetupSingletons(random_seed);
	SetupAndRunCartilageSheetSimulation(random_seed, random_birth_times, random_division_directions,
			cartesian_grid,
			n_cells_wide, n_cells_deep, n_cells_high, n_peri_lower, n_peri_upper,
			n_boundaries, upper_boundary, activation_percentage,
			maximum_perturbation, spring_stiffness, spring_stiffness_repulsion,
			homotypic_chondro_multiplier, baseline_adhesion_multiplier,
			simulation_end_time, force_function, output_directory);	

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

void SetForceFunction(OffLatticeSimulation<3>& simulator, std::string forceFunction,  
		double spring_stiffness, double spring_stiffness_repulsion, double alpha,
		double homotypic_peri_multiplier, double homotypic_chondro_multiplier, double heterotypic_multiplier){


	if (forceFunction.compare("cubic")==0){
		MAKE_PTR(CubicGeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		simulator.AddForce(p_force);
	}
	else if (forceFunction.compare("cubic_repulsion_only")==0){
		MAKE_PTR(RepulsionCubicForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		simulator.AddForce(p_force);
	}
	else if (forceFunction.compare("pwq")==0){
		MAKE_PTR(PWQGeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		p_force->SetRepulsionSpringStiffness(spring_stiffness_repulsion); // our default value fixed by experiments on optimal relative column height
		simulator.AddForce(p_force);
	}
	else if (forceFunction.compare("GLS_repulsion_only")==0){
		MAKE_PTR(RepulsionForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		simulator.AddForce(p_force);
	}
	else {
		// default is CellTissueTypeBasedGeneralisedSpringForce
		MAKE_PTR(CellTissueTypeBasedGeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		p_force->SetRepulsionSpringStiffness(spring_stiffness_repulsion); // our default value fixed by experiments on optimal relative column height
		p_force->SetAlpha(alpha);
		p_force->SetHomotypicPerichondrialSpringConstantMultiplier(
			homotypic_peri_multiplier);
		p_force->SetHomotypicChondrocyteSpringConstantMultiplier(
			homotypic_chondro_multiplier);
		p_force->SetHeterotypicSpringConstantMultiplier(heterotypic_multiplier);
		simulator.AddForce(p_force);
	}

}

void SetupAndRunCartilageSheetSimulation(unsigned random_seed,
		bool random_birth_times, bool random_division_directions, 
		bool cartesian_grid,
		unsigned n_cells_wide, unsigned n_cells_deep,
		unsigned n_cells_high, unsigned n_peri_lower, unsigned n_peri_upper,
		unsigned n_boundaries, double upper_boundary,
		double activation_percentage,
		double maximum_perturbation, double spring_stiffness,
		double spring_stiffness_repulsion,
		double homotypic_chondro_multiplier,
		double baseline_adhesion_multiplier, 
		double simulation_endtime,
		std::string force_function,
		std::string output_directory) {
	

	//bool random_birth_times = true;
	output_directory.append(boost::lexical_cast<std::string>(random_seed));

	double alpha;
	if (spring_stiffness == 0) {
		alpha = 1.0; // magic number, doesn't matter anyway
	} else {
		alpha = -2.0 * log(2.0 / spring_stiffness * 0.001); //not defined if spring_stiffness == 0
	}

	double homotypic_peri_multiplier = 1.0;
	//double homotypic_chondro_multiplier = 1.0;
	double heterotypic_multiplier = 0.5
			* (homotypic_peri_multiplier + homotypic_chondro_multiplier);

	CellBasedEventHandler::Enable();

	NodeBasedCartilageSheet* p_cartilage_sheet = new NodeBasedCartilageSheet();

	// set the sheet dimensions
	p_cartilage_sheet->SetCartilageSheetDimensions(n_cells_wide, n_cells_deep,
			n_cells_high);
	p_cartilage_sheet->setMaxCoordinatePerturbation(maximum_perturbation);

	if (!random_birth_times) {
		p_cartilage_sheet->setSynchronizeCellCycles(true);
	}
	if (random_division_directions){
		p_cartilage_sheet->setDivisionDirections(false);
	}

	if (cartesian_grid){
		p_cartilage_sheet->GenerateNodesOnCartesianGrid();
	}
	else {
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();
	}
	

	// setup the cell population
	if (!p_cartilage_sheet->isCellPopulationSetup()) {
		p_cartilage_sheet->Setup();
	}

	p_cartilage_sheet->setNumberOfPerichondrialLayersBelow(n_peri_lower);
	p_cartilage_sheet->setNumberOfPerichondrialLayersAbove(n_peri_upper);

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

	OffLatticeSimulation<3> simulator(*cell_population);
	//OffLatticeSimulation<3> simulator(cell_population);
	simulator.SetOutputDirectory(output_directory);
	simulator.SetEndTime(simulation_endtime); //hours
	simulator.SetSamplingTimestepMultiple(12);

	// call helper function to set force function
	SetForceFunction(simulator, force_function, 
						spring_stiffness, spring_stiffness_repulsion, 
						alpha, homotypic_peri_multiplier, 
						homotypic_chondro_multiplier, heterotypic_multiplier);

	if (n_boundaries > 0){
		//bottom plane
    	c_vector<double,3> point = zero_vector<double>(3);
    	c_vector<double,3> normal = zero_vector<double>(3);
    	normal(2) = -1.0;
		NodeBasedCellPopulation<3> nCellPop = *cell_population;
		//PlaneBoundaryCondition<3>* mypc = new PlaneBoundaryCondition<3>(&nCellPop, point, normal);
   		MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc, (cell_population.get(), point, normal));
    	//p_bc->SetUseJiggledNodesOnPlane(true);
    	simulator.AddCellPopulationBoundaryCondition(p_bc);

		if (n_boundaries == 2) {
			//upper plane
    		c_vector<double,3> point_up = zero_vector<double>(3);
    		point_up(2) = upper_boundary;
    		c_vector<double,3> normal_up = zero_vector<double>(3);
    		normal_up(2) = 1.0;
    		MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_up, (cell_population.get(), point_up, normal_up));
    		//p_bc->SetUseJiggledNodesOnPlane(true);
    		simulator.AddCellPopulationBoundaryCondition(p_bc_up);
		}
	}

	// Add the PatchSizeTracker to ensure that patches have maximum 6 cells
	MAKE_PTR(PatchSizeTrackingModifier<3>, p_modifier);
    simulator.AddSimulationModifier(p_modifier);
       

	CellBasedEventHandler::Reset();
	simulator.Solve();
	


	// write sheet parameters to file
	std::stringstream ss;
	ss << "/home/kubuntu1804/Documents/sf_simulation_results/" << output_directory << "/results_from_time_0/sheet.parameters";
	std::string sheet_params_filename = ss.str();
//	std::cout << sheet_params_filename << std::endl;
	std::ofstream sheet_params_file;
	sheet_params_file.open(sheet_params_filename.c_str());
//	std::cout << sheet_params_file.is_open() << std::endl;
	sheet_params_file << "---------------------------------\n";
	sheet_params_file << "Parameters of current sheet simulation:\n";
	sheet_params_file << "---------------------------------\n";
	sheet_params_file << "Random seed : " << random_seed << "\n";
	sheet_params_file << "Number cells X : " << n_cells_wide << "\n";
	sheet_params_file << "Number cells Y : " << n_cells_deep << "\n";
	sheet_params_file << "Number cells Z : " << n_cells_high << "\n";
	sheet_params_file << "Random birth times : " << random_birth_times << "\n";
	sheet_params_file << "Activation percentage : " << activation_percentage
			<< "\n";
	sheet_params_file << "Maximum perturbation : " << maximum_perturbation
			<< "\n";
//	sheet_params_file << "Adhesion spring stiffness : " << spring_stiffness << "\n";
//	sheet_params_file << "Attraction force decay : " << alpha << "\n";
	sheet_params_file << "Simulation end time : " << simulation_endtime << "\n";
//	sheet_params_file << "Force Function: " << force_function << "\n";
	sheet_params_file << "Output directory : " << output_directory << "\n";
	sheet_params_file.close();

//	std::cout << "Written sheet parameters to file." << std::endl;

	CellBasedEventHandler::Headings();
	CellBasedEventHandler::Report();
}



