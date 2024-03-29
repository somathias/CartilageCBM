/*
 * MesenchymalCondensationSimulation.cpp
 *
 *  Created on: Mar 12, 2020
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

#include "NodeBasedMesenchymalCondensation.hpp"
#include "OffLatticeSimulationDirectedDivision.hpp"
#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
#include "IndividualSpringStiffnessGeneralisedLinearSpringForce.hpp"
#include "CubicGeneralisedLinearSpringForce.hpp"
#include "RepulsionCubicForce.hpp"
#include "PWQGeneralisedLinearSpringForce.hpp"
#include "RepulsionForce.hpp"
#include "PatchSizeTrackingModifier.hpp"
#include "ChondrocytesOnlyCellCycleModel.hpp"


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
void SetupAndRunMesenchymalCondensationSimulation(unsigned randomSeed, bool, bool,  bool, bool, unsigned,
		unsigned, double, double, double, unsigned, double, double, double, double, double, double, double, double, std::string, std::string);
void SetForceFunction(OffLatticeSimulation<3>&, std::string, double, double, double, double, double, double);

int main(int argc, char *argv[]) {
	// This sets up PETSc and prints out copyright information, etc.
	//ExecutableSupport::StandardStartup(&argc, &argv);
	ExecutableSupport::InitializePetsc(&argc, &argv);

	// Define command line options
	boost::program_options::options_description general_options(
			"This is a Chaste executable.\n");
	general_options.add_options()("help", "Produce help message")("rdd",
			"Random division directions")("continue",
			"Continue simulation after increasing upper boundary and patch size limit")("flat",
			"Do not use offset in z direction")("lb0", "Have lower boundary at z=0")("S",
			boost::program_options::value<unsigned>()->default_value(0),
			"The random seed")("sw",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in x direction")("sd",
			boost::program_options::value<unsigned>()->default_value(5),
			"The number of cells in y direction")("sc",
			boost::program_options::value<double>()->default_value(1.0),
			"Scaling of the sheet in x-y direction")("u",
			boost::program_options::value<double>()->default_value(7.0),
			"The distance of the upper boundary to the lower one")("mu",
			boost::program_options::value<double>()->default_value(15.0),
			"The adhesion spring stiffness")("mu_R",
			boost::program_options::value<double>()->default_value(1.4),
			"The repulsion spring stiffness")("A",
			boost::program_options::value<double>()->default_value(0.5),
			"The percentage of activated stem cells")("psl",
			boost::program_options::value<unsigned>()->default_value(6),
			"The size limit for clonal patches")("p",
			boost::program_options::value<double>()->default_value(0.0),
			"The maximum perturbation of the initial coordinates.")("z", 
            boost::program_options::value<double>()->default_value(0.0), 
            "The maximum perturbation for the zenith angle if used." )("g1t",
			boost::program_options::value<double>()->default_value(30.0),
			"The transit cell G1 duration.")("ds",
			boost::program_options::value<double>()->default_value(10.0),
			"The s phase duration.")("T",
			boost::program_options::value<double>()->default_value(10.0),
			"The simulation end time")("dt",
			boost::program_options::value<double>()->default_value(0.008333),
			"The simulation time step")("F",
			boost::program_options::value<std::string>()->default_value(
					"pwq"), "The force function used")("output-dir",
			boost::program_options::value<std::string>()->default_value(
					"3dNodeBasedMesenchymalCondensation/"), "The output directory");

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

	// bool random_birth_times = true;
	// // set random birth times to false = synchronized case if wanted
	// if (variables_map.count("sbt")) {
	// 	random_birth_times = false;
	// }
	bool random_division_directions = false;
	if (variables_map.count("rdd")) {
		random_division_directions = true;
	}
	bool cont = false;
	if (variables_map.count("continue")) {
		cont = true;
	}

	bool use_offset = true;
	if (variables_map.count("flat")) {
		use_offset = false;
	}

	bool symmetrical_boundary = true;
	if (variables_map.count("lb0")) {
		symmetrical_boundary = false;
	}
	

	// Get ID and name from command line
	unsigned random_seed = variables_map["S"].as<unsigned>();
	unsigned n_cells_wide = variables_map["sw"].as<unsigned>();
	unsigned n_cells_deep = variables_map["sd"].as<unsigned>();
	double scaling = variables_map["sc"].as<double>();

	double upper_boundary = variables_map["u"].as<double>();
	double activation_percentage = variables_map["A"].as<double>();
	double maximum_perturbation = variables_map["p"].as<double>();
    double maximum_zenith_angle = variables_map["z"].as<double>();

	unsigned patch_size_limit = variables_map["psl"].as<unsigned>();

	double transit_cell_g1_duration = variables_map["g1t"].as<double>();
	double s_phase_duration = variables_map["ds"].as<double>();

	//std::cout <<"Patch size limit: " << patch_size_limit << std::endl;

	double spring_stiffness = variables_map["mu"].as<double>();
	double spring_stiffness_repulsion = variables_map["mu_R"].as<double>();
	double simulation_end_time = variables_map["T"].as<double>();
	double dt = variables_map["dt"].as<double>();

	std::string force_function =
			variables_map["F"].as<std::string>();
	std::string output_directory =
			variables_map["output-dir"].as<std::string>();


	SetupSingletons(random_seed);
	SetupAndRunMesenchymalCondensationSimulation(random_seed, random_division_directions, cont, use_offset,
			symmetrical_boundary,
			n_cells_wide, n_cells_deep, scaling, upper_boundary, activation_percentage, patch_size_limit,
			maximum_perturbation, maximum_zenith_angle, transit_cell_g1_duration, s_phase_duration, spring_stiffness, spring_stiffness_repulsion,
			simulation_end_time, dt, force_function, output_directory);	

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
		double spring_stiffness, double spring_stiffness_repulsion, double alpha){


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
        // if the string is not recognized just use GLS_repulsion_only
        MAKE_PTR(RepulsionForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		simulator.AddForce(p_force);
	}

}

void SetupAndRunMesenchymalCondensationSimulation(unsigned random_seed,
		bool random_division_directions, bool cont, bool use_offset, bool symmetrical_boundary,
		unsigned n_cells_wide, unsigned n_cells_deep, double scaling, double upper_boundary, 
		double activation_percentage, unsigned patch_size_limit,
		double maximum_perturbation, double maximum_zenith_angle, double transit_cell_g1_duration, double s_phase_duration,
		double spring_stiffness,
		double spring_stiffness_repulsion,
		double simulation_endtime, double dt,
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

	CellBasedEventHandler::Enable();

	NodeBasedMesenchymalCondensation* p_condensation = new NodeBasedMesenchymalCondensation();

	// set the sheet dimensions
	p_condensation->SetDimensions(n_cells_wide, n_cells_deep);
	p_condensation->setMaxCoordinatePerturbation(maximum_perturbation);
    p_condensation->setMaxZenithAnglePerturbation(maximum_zenith_angle);
	p_condensation->setDistanceBetweeenBoundaries(upper_boundary);
	p_condensation->SetPatchSizeLimit(patch_size_limit);	
	p_condensation->SetPhaseDurations(transit_cell_g1_duration, s_phase_duration);


	// if (!random_birth_times) {
	// 	p_condensation->setSynchronizeCellCycles(true);
	// }
	if (random_division_directions){
		p_condensation->setDivisionDirections(false);
	}
	// generate the nodes
	p_condensation->GenerateNodesOnHCPGrid(scaling, use_offset);

	// setup the cell population
	if (!p_condensation->isCellPopulationSetup()) {
		p_condensation->Setup();
	}

	// setup the initial stem cell configuration
	unsigned n_activated_stem_cells = floor(
			n_cells_wide * n_cells_deep * activation_percentage);
	p_condensation->InitialiseRandomConfiguration(
			n_activated_stem_cells);

	// get the cell population
	boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
			p_condensation->GetCellPopulation();

	OffLatticeSimulation<3> simulator(*cell_population);
	//OffLatticeSimulation<3> simulator(cell_population);
	simulator.SetOutputDirectory(output_directory);
	simulator.SetEndTime(simulation_endtime); //hours
	simulator.SetDt(dt);
	simulator.SetSamplingTimestepMultiple(12);

	// call helper function to set force function
	SetForceFunction(simulator, force_function, 
						spring_stiffness, spring_stiffness_repulsion, 
						alpha);

	//bottom plane
    c_vector<double,3> point = zero_vector<double>(3);
	if(symmetrical_boundary){
		point(2) = - upper_boundary/2.0;
	}
	c_vector<double,3> normal = zero_vector<double>(3);
    normal(2) = -1.0;
	//NodeBasedCellPopulation<3> nCellPop = *cell_population;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_down, (cell_population.get(), point, normal));
    //p_bc->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc_down);

    //upper plane
    c_vector<double,3> point_up = zero_vector<double>(3);
	if(symmetrical_boundary){
		point_up(2) = upper_boundary/2.0;
	}
	else {
    	point_up(2) = upper_boundary;
	}
    c_vector<double,3> normal_up = zero_vector<double>(3);
    normal_up(2) = 1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_up, (cell_population.get(), point_up, normal_up));
    //p_bc->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc_up);

	// Add the PatchSizeTracker to ensure that patches have maximum 6 cells
	MAKE_PTR(PatchSizeTrackingModifier<3>, p_modifier);
    simulator.AddSimulationModifier(p_modifier);
   

	CellBasedEventHandler::Reset();
	simulator.Solve();

	if (cont){
		//reset the both the upper and the lower plane
		point(2) = -1.0;
		point_up(2) = upper_boundary + 1.0;
		simulator.RemoveAllCellPopulationBoundaryConditions();
		MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_lower, (cell_population.get(), point, normal));
    	simulator.AddCellPopulationBoundaryCondition(p_bc_lower);
		MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_up_wider, (cell_population.get(), point_up, normal_up));
    	simulator.AddCellPopulationBoundaryCondition(p_bc_up_wider);

		//increase the patch size limit
		for (AbstractCellPopulation<3>::Iterator cell_iter =
			cell_population.get()->Begin(); cell_iter != cell_population.get()->End();
			++cell_iter) {
			
			static_cast<ChondrocytesOnlyCellCycleModel*>(cell_iter->GetCellCycleModel())->SetPatchSizeLimit(patch_size_limit+2);
			cell_iter->SetBirthTime(simulation_endtime); //reset the birth time so that not every cell wants to divide at once
		}
		simulator.SetEndTime(simulation_endtime+50); //hours

		simulator.Solve();
	}

	

	


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
	sheet_params_file << "Number cells Z : " << 1 << "\n";
	// sheet_params_file << "Random birth times : " << random_birth_times << "\n";
	sheet_params_file << "Activation percentage : " << activation_percentage
			<< "\n";
	sheet_params_file << "Maximum perturbation : " << maximum_perturbation
			<< "\n";
    sheet_params_file << "Maximum zenith angle : " << maximum_zenith_angle
			<< "\n";
	sheet_params_file << "Distance between upper and lower boundary : " << upper_boundary
			<< "\n";
	sheet_params_file << "Patch Size Limit : " << patch_size_limit  << "\n";
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



