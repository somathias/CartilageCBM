#include <cxxtest/TestSuite.h>
#include <iostream>
#include <sstream>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "NodeBasedMesenchymalCondensation.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"

#include "PatchSizeTrackingModifier.hpp"

#include "FakePetscSetup.hpp"

class TestNodeBasedMesenchymalCondensation: public AbstractCellBasedTestSuite {
public:
	void TestSheet()  {

        EXIT_IF_PARALLEL;

		// Construct a new mesenchymal condensation
		NodeBasedMesenchymalCondensation* p_condensation =
				new NodeBasedMesenchymalCondensation();

		// set the sheet dimensions
		p_condensation->SetDimensions(5, 4);
		p_condensation->setMaxCoordinatePerturbation(0.1);
		p_condensation->UseRandomSeed();
		unsigned seed = p_condensation->getSeed();
		std::stringstream ss;
		ss << "/" << seed;
		std::string seed_string = ss.str();

        std::cout <<"Before generating the nodes"<< std::endl;
		// generate the nodes
		p_condensation->GenerateNodesOnCartesianGrid();

        std::cout <<"After generating the nodes"<< std::endl;


		// setup the cell population
		if (!p_condensation->isCellPopulationSetup()) {
			p_condensation->Setup();
		}

        std::cout <<"After setting up the population"<< std::endl;


		// setup the initial transit cell configuration
		//p_condensation->InitialiseBulktransitCellConfiguration(2, 1);
		p_condensation->InitialiseRandomConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
				p_condensation->GetCellPopulation();

		//pass it to the simulator
		OffLatticeSimulation<3> simulator(*cell_population);
		simulator.SetOutputDirectory("NodeBasedMesenchymalCondensation" + seed_string);
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(50.0);

		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(1.0);
		simulator.AddForce(p_force);

		std::cout <<"After adding force"<< std::endl;

		// Add the PatchSizeTracker to ensure that patches have maximum 6 cells
		MAKE_PTR(PatchSizeTrackingModifier<3>, p_modifier);
    	simulator.AddSimulationModifier(p_modifier);

		std::cout <<"After adding modifier"<< std::endl;



		simulator.Solve();
		std::cout <<"After solve"<< std::endl;

	}

	void TestRandomTransitCellConfiguration()  {

		EXIT_IF_PARALLEL;

		// Construct a new cartilage sheet
		NodeBasedMesenchymalCondensation* p_condensation =
		new NodeBasedMesenchymalCondensation();

		// set the sheet dimensions
		p_condensation->SetDimensions(5, 4);
		// generate the nodes
		p_condensation->GenerateNodesOnCartesianGrid();

		// setup the cell population
		p_condensation->Setup();

		// setup the initial activated cell configuration
		p_condensation->InitialiseRandomConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
		p_condensation->GetCellPopulation();


		unsigned n_transit_cells = 0;
		unsigned n_diff_cells = 0;
		for (AbstractCellPopulation<3>::Iterator cell_iter =
					cell_population->Begin(); cell_iter != cell_population->End();
					++cell_iter) {
			if (cell_iter->GetCellProliferativeType()->IsType<TransitCellProliferativeType>()) {
				n_transit_cells++;
			}
			else if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()){
				n_diff_cells++;
			}
			

		}
		TS_ASSERT_EQUALS(n_transit_cells, 5);//number of transit cells should be 5;
		TS_ASSERT_EQUALS(n_diff_cells, 15);//number of differentiated cells should be 5*4-5;

	}

	void TestPatchSizeTrackingModifier(){
		EXIT_IF_PARALLEL;

		// Construct a new mesenchymal condensation
		NodeBasedMesenchymalCondensation* p_condensation =
				new NodeBasedMesenchymalCondensation();

		// set the sheet dimensions
		p_condensation->SetDimensions(5, 4);
		p_condensation->setMaxCoordinatePerturbation(0.1);
		p_condensation->UseRandomSeed();
		unsigned seed = p_condensation->getSeed();
		std::stringstream ss;
		ss << "/" << seed;
		std::string seed_string = ss.str();

        std::cout <<"Before generating the nodes"<< std::endl;
		// generate the nodes
		p_condensation->GenerateNodesOnCartesianGrid();

        std::cout <<"After generating the nodes"<< std::endl;


		// setup the cell population
		if (!p_condensation->isCellPopulationSetup()) {
			p_condensation->Setup();
		}

        std::cout <<"After setting up the population"<< std::endl;


		// setup the initial transit cell configuration
		//p_condensation->InitialiseBulktransitCellConfiguration(2, 1);
		p_condensation->InitialiseRandomConfiguration(5);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
				p_condensation->GetCellPopulation();

		//pass it to the simulator
		OffLatticeSimulation<3> simulator(*cell_population);
		simulator.SetOutputDirectory("NodeBasedMesenchymalCondensation" + seed_string);
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(150.0);

		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(1.0);
		simulator.AddForce(p_force);

	
		// Add the PatchSizeTracker to ensure that patches have maximum 6 cells
		MAKE_PTR(PatchSizeTrackingModifier<3>, p_modifier);
    	simulator.AddSimulationModifier(p_modifier);

		simulator.Solve();
	}

	/**
	 * Minimal testing for the generation of the node coordinates on a cartesian lattice
	 */
	void TestCartesianNodeGeneration()  {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 2;

		// Construct a new cartilage sheet
		NodeBasedMesenchymalCondensation* p_condensation =
		new NodeBasedMesenchymalCondensation();

		// set the sheet dimensions
		p_condensation->SetDimensions(n_nodes_width,
				n_nodes_depth);
		//p_condensation->setMaxCoordinatePerturbation(0.1);
		//p_condensation->UseRandomSeed();
		// generate the nodes
		p_condensation->GenerateNodesOnCartesianGrid();

		TS_ASSERT_EQUALS(p_condensation->mNodes.size(),
				n_nodes_width * n_nodes_depth);

		for (unsigned i = 0; i < p_condensation->mNodes.size(); i++) {
			c_vector<double, 3> coordinates =
			p_condensation->mNodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.0);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.0);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[2], 7.0); //Right now the z-value is drawn uniformly from [0,7]
		}

		// now with perturbation
		p_condensation->setMaxCoordinatePerturbation(0.5);
		p_condensation->GenerateNodesOnCartesianGrid();

		for (unsigned i = 0; i < p_condensation->mNodes.size(); i++) {
			c_vector<double, 3> coordinates =
			p_condensation->mNodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[2], 7.5); //Right now the z-value is drawn uniformly from [0,7]
		}

	}

    /**
	 * Minimal testing for the generation of the node coordinates on a HCP lattice
	 */
	void TestHCPNodeGeneration()  {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 3;

		// Construct a new cartilage sheet
		NodeBasedMesenchymalCondensation* p_condensation =
		new NodeBasedMesenchymalCondensation();

		// set the sheet dimensions
		p_condensation->SetDimensions(n_nodes_width,
				n_nodes_depth);
		//p_condensation->setMaxCoordinatePerturbation(0.1);
		//p_condensation->UseRandomSeed();
		// generate the nodes
		p_condensation->GenerateNodesOnHCPGrid();

		TS_ASSERT_EQUALS(p_condensation->mNodes.size(),
				n_nodes_width * n_nodes_depth);

		c_vector<double, 3> coordinates_first =
		p_condensation->mNodes[0]->rGetLocation();
		c_vector<double, 3> coordinates_second =
		p_condensation->mNodes[1]->rGetLocation();
		c_vector<double, 3> coordinates_last =
		p_condensation->mNodes[p_condensation->mNodes.size() - 1]->rGetLocation();

		//check that first node is in origin
		TS_ASSERT_EQUALS(coordinates_first[0], 0);
		TS_ASSERT_EQUALS(coordinates_first[1], 0);
		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[2], 7.0);

		// //check that first two have distance 1 cell diameter
		// double distance =
		// sqrt(
		// 		(coordinates_first[0] - coordinates_second[0])
		// 		* (coordinates_first[0] - coordinates_second[0])
		// 		+ (coordinates_first[1] - coordinates_second[1])
		// 		* (coordinates_first[1]
		// 				- coordinates_second[1])
		// 		+ (coordinates_first[2] - coordinates_second[2])
		// 		* (coordinates_first[2]
		// 				- coordinates_second[2]));
		// TS_ASSERT_DELTA(distance, 1, 1e-4);

		// //check the coordinates of the last node
		// TS_ASSERT_EQUALS(coordinates_last[0], 2);
		// TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-4);
		// TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-4);

		// now check with perturbation
		p_condensation->setMaxCoordinatePerturbation(0.1);
		p_condensation->GenerateNodesOnHCPGrid();

		// coordinates_first = p_condensation->mNodes[0]->rGetLocation();
		// coordinates_second = p_condensation->mNodes[1]->rGetLocation();
		// coordinates_last =
		// p_condensation->mNodes[p_condensation->mNodes.size() - 1]->rGetLocation();

		// //check that first node is in origin (+perturbation)
		// TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[0], 0.1);
		// TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[1], 0.1);
		// TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[2], 7.1);

		// //check that first two have distance 1 cell diameter
		// distance =
		// sqrt(
		// 		(coordinates_first[0] - coordinates_second[0])
		// 		* (coordinates_first[0] - coordinates_second[0])
		// 		+ (coordinates_first[1] - coordinates_second[1])
		// 		* (coordinates_first[1]
		// 				- coordinates_second[1])
		// 		+ (coordinates_first[2] - coordinates_second[2])
		// 		* (coordinates_first[2]
		// 				- coordinates_second[2]));
		// TS_ASSERT_LESS_THAN_EQUALS(abs(distance - 1), 0.2);

		// //check the coordinates of the last node
		// TS_ASSERT_LESS_THAN_EQUALS(coordinates_last[0], 2.1);
		// TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-1);
		// TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-1);
	}

};
