/*
 * TestIndividualSpringStiffnessGeneralisedLinearSpringForce.hpp
 *
 *  Created on: Jan 18, 2018
 *      Author: Sonja Mathias
 */

#ifndef TESTINDIVIDUALSPRINGSTIFFNESSGENERALISEDLINEARSPRINGFORCE_HPP_
#define TESTINDIVIDUALSPRINGSTIFFNESSGENERALISEDLINEARSPRINGFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "IndividualSpringStiffnessGeneralisedLinearSpringForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestIndividualSpringStiffnessGeneralisedLinearSpringForce: public AbstractCellBasedTestSuite {
public:

	void TestIndividualSpringStiffnessGeneralisedLinearSpringForceWithMeshBasedCellPopulation()
			 {
		EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

		unsigned cells_across = 7;
		unsigned cells_up = 5;
		unsigned thickness_of_ghost_layer = 3;

		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0,1);

		HoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		// Create cells
		std::vector<CellPtr> cells;
		CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);

		// Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		// Create force
		IndividualSpringStiffnessGeneralisedLinearSpringForce<2> force;

		// Test get method
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(), 15.0, 1e-6);
		TS_ASSERT_DELTA(force.GetAlpha(), 5.0, 1e-6);


		// Initialise a vector of node forces
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the x-axis and calculate the force exerted on a neighbour
		c_vector<double,2> old_point = p_mesh->GetNode(59)->rGetLocation();
		ChastePoint<2> new_point;
		new_point.rGetLocation()[0] = old_point[0]+0.5;
		new_point.rGetLocation()[1] = old_point[1];
		p_mesh->SetNode(59, new_point, false);

		double spring_stiffness = force.GetMeinekeSpringStiffness();

		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}
		force.AddForceContribution(cell_population);

		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], (-3+4.0/sqrt(7.0))*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

	}

	void TestIndividualSpringStiffnessGeneralisedLinearSpringForceWithNodeBasedCellPopulation() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, -1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u, false, 1.0, 0.0, 0.0));
//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		// Create force
		IndividualSpringStiffnessGeneralisedLinearSpringForce<3> force;

		// Test set/get method
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(), 15.0, 1e-6);
		TS_ASSERT_DELTA(force.GetAlpha(), 5.0, 1e-6);

		double alpha = 12.5;
		force.SetAlpha(alpha);
		TS_ASSERT_DELTA(force.GetAlpha(), alpha, 1e-6);

		double repulsion_spring_stiffness = 1.0;
		force.SetRepulsionSpringStiffness(repulsion_spring_stiffness);
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(), repulsion_spring_stiffness, 1e-6);

		// Initialise a vector of node forces
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the x-axis and calculate the force exerted on a neighbour
		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_old_coords = cell_population.GetLocationOfCellCentre(my_cell);
		ChastePoint<3> new_coords;
		new_coords.rGetLocation()[0] = my_old_coords[0]+0.25;
		new_coords.rGetLocation()[1] = my_old_coords[1];
		new_coords.rGetLocation()[2] = my_old_coords[2];
		unsigned node_index = cell_population.GetLocationIndexUsingCell(my_cell);
		cell_population.SetNode(node_index, new_coords);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], 0.25, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[2], 0.0, 1e-4);

		double spring_stiffness_adhesion = force.GetMeinekeSpringStiffness();
		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution = force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], spring_stiffness_repulsion*log(0.75), 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], -0.25*spring_stiffness_adhesion*exp(-alpha*0.25), 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		cell_population.Update(false);
		std::vector< std::pair<Node<3>*, Node<3>* > >& r_node_pairs = cell_population.rGetNodePairs();
		TS_ASSERT_EQUALS(r_node_pairs.size(), 3);

		// Initialise a vector of node forces
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}
		force.AddForceContribution(cell_population);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], spring_stiffness_repulsion*log(0.75) - 0.25*spring_stiffness_adhesion*exp(-alpha*0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], 0.25*spring_stiffness_adhesion*exp(-alpha*0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], -spring_stiffness_repulsion*log(0.75), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0, 1e-4);

	}
	void TestIndividualSpringStiffnessGeneralisedLinearSpringForceOutputParameters()
	{
		EXIT_IF_PARALLEL;
		std::string output_directory = "TestForcesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with IndividualSpringStiffnessGeneralisedLinearSpringForce
		IndividualSpringStiffnessGeneralisedLinearSpringForce<2> linear_force;
		linear_force.SetCutOffLength(1.5);
		TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "IndividualSpringStiffnessGeneralisedLinearSpringForce-2-2");

		out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("individual_spring_stiffness_results.parameters");
		linear_force.OutputForceParameters(linear_force_parameter_file);
		linear_force_parameter_file->close();

		{
			FileFinder generated_file = output_file_handler.FindFile("individual_spring_stiffness_results.parameters");
			FileFinder reference_file("projects/cartilage/test/data/TestIndividualSpringStiffnessGeneralisedLinearSpringForce/individual_spring_stiffness_results.parameters",
					RelativeTo::ChasteSourceRoot);
			FileComparison comparer(generated_file,reference_file);
			TS_ASSERT(comparer.CompareFiles());
		}

	}

	void TestIndividualSpringStiffnessGeneralisedLinearSpringForceArchiving() 
	{
		EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "IndividualSpringStiffnessGeneralisedLinearSpringForce.arch";

		{
			IndividualSpringStiffnessGeneralisedLinearSpringForce<2> force;

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Set member variables
			force.SetMeinekeSpringStiffness(12.34);
			force.SetMeinekeDivisionRestingSpringLength(0.856);
			force.SetMeinekeSpringGrowthDuration(2.593);
			force.SetRepulsionSpringStiffness(15.67);
			force.SetAlpha(12.5);

			// Serialize via pointer to most abstract class possible
			AbstractForce<2>* const p_force = &force;
			output_arch << p_force;
		}

		{
			AbstractForce<2>* p_force;

			// Create an input archive
			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			// Restore from the archive
			input_arch >> p_force;

			// Test member variables
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetRepulsionSpringStiffness(), 15.67, 1e-6);
			TS_ASSERT_DELTA((static_cast<IndividualSpringStiffnessGeneralisedLinearSpringForce<2>*>(p_force))->GetAlpha(), 12.5, 1e-6);

			// Tidy up
			delete p_force;
		}
	}

};



#endif /* TESTINDIVIDUALSPRINGSTIFFNESSGENERALISEDLINEARSPRINGFORCE_HPP_ */
