#include <cxxtest/TestSuite.h>
#include <iostream>
#include <sstream>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "NodeBasedCartilageSheet.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "FakePetscSetup.hpp"

class TestNodeBasedCartilageSheet: public AbstractCellBasedTestSuite {
public:
	void TestSheet()  {

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
				new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(5, 4, 3);
		p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		p_cartilage_sheet->UseRandomSeed();
		unsigned seed = p_cartilage_sheet->getSeed();
		std::stringstream ss;
		ss << "/" << seed;
		std::string seed_string = ss.str();

		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		// setup the cell population
		if (!p_cartilage_sheet->isCellPopulationSetup()) {
			p_cartilage_sheet->Setup();
		}
		// setup the cell tissue types and cell division directions
		p_cartilage_sheet->InitialiseTissueLayersAndCellDivisionDirections();
		// setup the initial stem cell configuration
		//p_cartilage_sheet->InitialiseBulkStemCellConfiguration(2, 1);
		p_cartilage_sheet->InitialiseRandomStemCellConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
				p_cartilage_sheet->GetCellPopulation();

		//pass it to the simulator
		OffLatticeSimulation<3> simulator(*cell_population);
		simulator.SetOutputDirectory("NodeBasedCartilageSheet" + seed_string);
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(50.0);

		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(1.0);
		simulator.AddForce(p_force);

		simulator.Solve();
	}

	void TestRandomStemCellConfiguration()  {

		EXIT_IF_PARALLEL;

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
		new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(5, 4, 2);
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		// setup the cell population
		p_cartilage_sheet->Setup();

		// setup the cell tissue types and cell division directions
		//p_cartilage_sheet->InitialiseTissueLayersAndCellDivisionDirections();
		// setup the initial stem cell configuration
		p_cartilage_sheet->InitialiseRandomStemCellConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
		p_cartilage_sheet->GetCellPopulation();


		unsigned n_stem_cells = 0;
		unsigned n_diff_cells = 0;
		for (AbstractCellPopulation<3>::Iterator cell_iter =
					cell_population->Begin(); cell_iter != cell_population->End();
					++cell_iter) {
			if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) {
				n_stem_cells++;
			}
			else if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()){
				n_diff_cells++;
			}

		}
		TS_ASSERT_EQUALS(n_stem_cells, 5);//number of stem cells should be 5;
		TS_ASSERT_EQUALS(n_diff_cells, 35);//number of differentiated cells should be 5*4*2-5;

	}

	void TestRandomStemCellConfigurationSingleLayer()  {

		EXIT_IF_PARALLEL;

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
		new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(5, 4, 1);
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		// setup the cell population
		p_cartilage_sheet->Setup();

		// setup the cell tissue types and cell division directions
		//p_cartilage_sheet->InitialiseTissueLayersAndCellDivisionDirections();
		// setup the initial stem cell configuration
		p_cartilage_sheet->InitialiseRandomStemCellConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
		p_cartilage_sheet->GetCellPopulation();


		unsigned n_stem_cells = 0;
		unsigned n_diff_cells = 0;
		for (AbstractCellPopulation<3>::Iterator cell_iter =
					cell_population->Begin(); cell_iter != cell_population->End();
					++cell_iter) {
			if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>()) {
				n_stem_cells++;
			}
			else if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()){
				n_diff_cells++;
			}

		}
		TS_ASSERT_EQUALS(n_stem_cells, 5);//number of stem cells should be 5;
		TS_ASSERT_EQUALS(n_diff_cells, 15);//number of differentiated cells should be 5*4*1-5;

	}

	/**
	 * Minimal testing for the generation of the node coordinates on a cartesian lattice
	 */
	void TestCartesianNodeGeneration()  {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 2;
		unsigned n_nodes_height = 1;

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
		new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(n_nodes_width,
				n_nodes_depth, n_nodes_height);
		//p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		//p_cartilage_sheet->UseRandomSeed();
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnCartesianGrid();

		TS_ASSERT_EQUALS(p_cartilage_sheet->mNodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		for (unsigned i = 0; i < p_cartilage_sheet->mNodes.size(); i++) {
			c_vector<double, 3> coordinates =
			p_cartilage_sheet->mNodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.0);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.0);
			TS_ASSERT_EQUALS(coordinates[2], 0.0);
		}

		// now with perturbation
		p_cartilage_sheet->setMaxCoordinatePerturbation(0.5);
		p_cartilage_sheet->GenerateNodesOnCartesianGrid();

		for (unsigned i = 0; i < p_cartilage_sheet->mNodes.size(); i++) {
			c_vector<double, 3> coordinates =
			p_cartilage_sheet->mNodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[2], 0.5);
		}

	}

	/**
	 * Minimal testing for the generation of the node coordinates on a HCP lattice
	 */
	void TestHCPNodeGeneration()  {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 3;
		unsigned n_nodes_height = 3;

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
		new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(n_nodes_width,
				n_nodes_depth, n_nodes_height);
		//p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		//p_cartilage_sheet->UseRandomSeed();
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		TS_ASSERT_EQUALS(p_cartilage_sheet->mNodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		c_vector<double, 3> coordinates_first =
		p_cartilage_sheet->mNodes[0]->rGetLocation();
		c_vector<double, 3> coordinates_second =
		p_cartilage_sheet->mNodes[1]->rGetLocation();
		c_vector<double, 3> coordinates_last =
		p_cartilage_sheet->mNodes[p_cartilage_sheet->mNodes.size() - 1]->rGetLocation();

		//check that first node is in origin
		TS_ASSERT_EQUALS(coordinates_first[0], 0);
		TS_ASSERT_EQUALS(coordinates_first[1], 0);
		TS_ASSERT_EQUALS(coordinates_first[2], 0);

		//check that first two have distance 1 cell diameter
		double distance =
		sqrt(
				(coordinates_first[0] - coordinates_second[0])
				* (coordinates_first[0] - coordinates_second[0])
				+ (coordinates_first[1] - coordinates_second[1])
				* (coordinates_first[1]
						- coordinates_second[1])
				+ (coordinates_first[2] - coordinates_second[2])
				* (coordinates_first[2]
						- coordinates_second[2]));
		TS_ASSERT_DELTA(distance, 1, 1e-4);

		//check the coordinates of the last node
		TS_ASSERT_EQUALS(coordinates_last[0], 2);
		TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-4);
		TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-4);

		// now check with perturbation
		p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		coordinates_first = p_cartilage_sheet->mNodes[0]->rGetLocation();
		coordinates_second = p_cartilage_sheet->mNodes[1]->rGetLocation();
		coordinates_last =
		p_cartilage_sheet->mNodes[p_cartilage_sheet->mNodes.size() - 1]->rGetLocation();

		//check that first node is in origin (+perturbation)
		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[0], 0.1);
		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[1], 0.1);
		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[2], 0.1);

		//check that first two have distance 1 cell diameter
		distance =
		sqrt(
				(coordinates_first[0] - coordinates_second[0])
				* (coordinates_first[0] - coordinates_second[0])
				+ (coordinates_first[1] - coordinates_second[1])
				* (coordinates_first[1]
						- coordinates_second[1])
				+ (coordinates_first[2] - coordinates_second[2])
				* (coordinates_first[2]
						- coordinates_second[2]));
		TS_ASSERT_LESS_THAN_EQUALS(abs(distance - 1), 0.2);

		//check the coordinates of the last node
		TS_ASSERT_LESS_THAN_EQUALS(coordinates_last[0], 2.1);
		TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-1);
		TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-1);
	}

};
