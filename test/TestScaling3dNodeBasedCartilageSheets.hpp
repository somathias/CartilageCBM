#include <ctime>
#include <string>
#include <sstream>
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
//#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
//#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
//#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellAncestorWriter.hpp"

#include "OffLatticeSimulationDirectedDivision.hpp"
#include "CellTissueTypeBasedCellCycleModel.hpp"
#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"
#include "NodeBasedCartilageSheet.hpp"

#include "PetscSetupAndFinalize.hpp"
/**
 * Second try to build a model for the cartilage sheet based on a center-based model.
 */

class TestScaling3dNodeBasedCartilageSheets: public AbstractCellBasedTestSuite {
public:

	/**
	 * Tests the 3D node-based cartilage sheet simulation
	 */
	void Test3dNodeBasedCartilageSheet() throw (Exception) {

		bool random_seed = true;
		unsigned n_cells_wide = 3;
		unsigned n_cells_deep = 3;
		unsigned n_cells_high = 3;
		unsigned n_differentiated_cells_width = 0;
		unsigned n_differentiated_cells_depth = 0;
		bool random_birth_times = true;

		double spring_stiffness = 1.0;
		double perichondrial_spring_constant_multiplier = 1.0;
		double chondrocyte_spring_constant_multiplier = 1.0;
		double heterotypic_spring_constant_multiplier = 1.0;
		double simulation_endtime = 20.0;

		std::string output_directory =
				"3dNodeBasedCartilageSheet/Test3dHCPCoordinates/";

		RunNodeBasedCartilageSheet(random_seed, n_cells_wide, n_cells_deep,
				n_cells_high, n_differentiated_cells_width,
				n_differentiated_cells_depth, random_birth_times,
				spring_stiffness, perichondrial_spring_constant_multiplier,
				chondrocyte_spring_constant_multiplier,
				heterotypic_spring_constant_multiplier, output_directory,
				simulation_endtime);

	}





private:

	/**
	 * Sets up the simulation for the specified input parameters
	 */
	void RunNodeBasedCartilageSheet(bool random_seed, unsigned n_cells_wide,
			unsigned n_cells_deep, unsigned n_cells_high,
			unsigned n_differentiated_cells_width,
			unsigned n_differentiated_cells_depth, bool random_birth_times,
			double spring_stiffness,
			double perichondrial_spring_constant_multiplier,
			double chondrocyte_spring_constant_multiplier,
			double heterotypic_spring_constant_multiplier,
			std::string output_directory, double simulation_endtime)
					throw (Exception) {
		CellBasedEventHandler::Enable();

		NodeBasedCartilageSheet* p_cartilage_sheet =
				new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(n_cells_wide,
				n_cells_deep, n_cells_high);
		p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		if (random_seed) {
			p_cartilage_sheet->UseRandomSeed();
			unsigned seed = p_cartilage_sheet->getSeed();
			std::stringstream ss;
			ss << "/" << seed;
			std::string seed_string = ss.str();
			output_directory.append(seed_string);
		}
		if (!random_birth_times){
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
		p_cartilage_sheet->InitialiseBulkStemCellConfiguration(
				n_cells_wide - 2 * n_differentiated_cells_width,
				n_cells_deep - 2 * n_differentiated_cells_depth);

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
		p_force->SetHomotypicPerichondrialSpringConstantMultiplier(
				perichondrial_spring_constant_multiplier);
		p_force->SetHomotypicChondrocyteSpringConstantMultiplier(
				chondrocyte_spring_constant_multiplier);
		p_force->SetHeterotypicSpringConstantMultiplier(
				heterotypic_spring_constant_multiplier);
		simulator.AddForce(p_force);

		CellBasedEventHandler::Reset();
		simulator.Solve();

		CellBasedEventHandler::Headings();
		CellBasedEventHandler::Report();

	}



};
