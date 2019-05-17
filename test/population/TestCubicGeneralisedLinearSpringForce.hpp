/*
 * TestCubicGeneralisedLinearSpringForce.hpp
 *
 *  Created on: Apr 30, 2019
 *      Author: Sonja Mathias
 */

#ifndef TESTCUBICGENERALISEDLINEARSPRINGFORCE_HPP_
#define TESTCUBICGENERALISEDLINEARSPRINGFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CubicGeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCubicGeneralisedLinearSpringForce: public AbstractCellBasedTestSuite {
public:

	void TestCubicGeneralisedLinearSpringForceWithMeshBasedCellPopulation()
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
		CubicGeneralisedLinearSpringForce<2> force;


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

	void TestCubicGeneralisedLinearSpringForceWithNodeBasedCellPopulation() {

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
		CubicGeneralisedLinearSpringForce<3> force;
        force.SetCutOffLength(1.5);

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

		double spring_stiffness = force.GetMeinekeSpringStiffness();

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution = force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], -spring_stiffness*9/64.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], -spring_stiffness/64.0, 1e-4);
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

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], -spring_stiffness*9/64.0 -spring_stiffness/64.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], spring_stiffness/64.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], spring_stiffness*9/64.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0, 1e-4);

	}



};



#endif /* TESTCUBICGENERALISEDLINEARSPRINGFORCE_HPP_ */
