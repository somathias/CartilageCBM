/*
 * TestDirectionalAdhesionForce.hpp
 *
 *  Created on: Mar 1, 2018
 *      Author: Sonja Mathias
 */

#ifndef TESTDIRECTIONALADHESIONFORCE_HPP_
#define TESTDIRECTIONALADHESIONFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "DirectionalAdhesionForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDirectionalAdhesionForce: public AbstractCellBasedTestSuite {
public:

	void TestDirectionalAdhesionForceWithNodeBasedCellPopulation() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, -1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u, false, 1.0, 0.0, 0.0));
		//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),
				p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		boost::shared_ptr<AbstractCellProperty> p_perichondrial(
				cell_population.GetCellPropertyRegistry()->Get<
						PerichondrialCellTissueType>());
		boost::shared_ptr<AbstractCellProperty> p_chondrocyte(
				cell_population.GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>());

		cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(
				p_perichondrial); //the cell in the middle is perichondrial
		cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(
				p_chondrocyte); // the other two are chondrocytes
		cell_population.GetCellUsingLocationIndex(2)->AddCellProperty(
				p_chondrocyte);

		// Create force
		DirectionalAdhesionForce force;

		// Test default values
		TS_ASSERT_DELTA(force.GetBaselineAdhesionMultiplier(), 0.1, 1e-6);

		// Test set/get method
		force.SetBaselineAdhesionMultiplier(0.0);

		TS_ASSERT_DELTA(force.GetBaselineAdhesionMultiplier(), 0.0, 1e-6);

		force.SetBaselineAdhesionMultiplier(0.1);
		double heterotypic_multiplier = 4.0;
		force.SetHeterotypicSpringConstantMultiplier(heterotypic_multiplier);

		double alpha = 12.5;
		force.SetAlpha(alpha);
		TS_ASSERT_DELTA(force.GetAlpha(), alpha, 1e-6);

		double repulsion_spring_stiffness = 1.0;
		force.SetRepulsionSpringStiffness(repulsion_spring_stiffness);
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(),
				repulsion_spring_stiffness, 1e-6);

		// Initialise a vector of node forces
		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the x-axis and calculate the force exerted on a neighbour
		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_old_coords =
				cell_population.GetLocationOfCellCentre(my_cell);
		ChastePoint<3> new_coords;
		new_coords.rGetLocation()[0] = my_old_coords[0] + 0.25;
		new_coords.rGetLocation()[1] = my_old_coords[1];
		new_coords.rGetLocation()[2] = my_old_coords[2];
		unsigned node_index = cell_population.GetLocationIndexUsingCell(
				my_cell);
		cell_population.SetNode(node_index, new_coords);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], 0.25,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[2], 0.0,
				1e-4);

		double spring_stiffness_adhesion = force.GetMeinekeSpringStiffness();
		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution =
				force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				spring_stiffness_repulsion * log(0.75), 1e-4); //repulsion force, hence no multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				-0.25 * heterotypic_multiplier * spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4); //adhesion force, hence heterotypic multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		cell_population.Update(false);
		std::vector<std::pair<Node<3>*, Node<3>*> >& r_node_pairs =
				cell_population.rGetNodePairs();
		TS_ASSERT_EQUALS(r_node_pairs.size(), 3);

		// Initialise a vector of node forces
		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
			cell_population.GetNode(i)->ClearAppliedForce();
		}
		force.AddForceContribution(cell_population);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0],
				spring_stiffness_repulsion * log(0.75)
						- 0.25 * heterotypic_multiplier
								* spring_stiffness_adhesion
								* exp(-alpha * 0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0],
				0.25 * heterotypic_multiplier * spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0],
				-spring_stiffness_repulsion * log(0.75), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0,
				1e-4);

	}

	void TestDirectionalAdhesionForcePerichondrialCase() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, -1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u, false, 1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(3u, false, 0.25, 0.0, 1.25));
		nodes.push_back(new Node<3>(4u, false, 1.25, 0.0, 1.0));
		//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),
				p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		boost::shared_ptr<AbstractCellProperty> p_perichondrial(
				cell_population.GetCellPropertyRegistry()->Get<
						PerichondrialCellTissueType>());

		// all cells are perichondrial cells
		cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(
				p_perichondrial);
		cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(
				p_perichondrial);
		cell_population.GetCellUsingLocationIndex(2)->AddCellProperty(
				p_perichondrial);
		cell_population.GetCellUsingLocationIndex(3)->AddCellProperty(
				p_perichondrial);
		cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(
				p_perichondrial);

		// Create force
		DirectionalAdhesionForce force;

		// Test default values
		double baseline_adhesion_multiplier = force.GetBaselineAdhesionMultiplier();
		TS_ASSERT_DELTA(baseline_adhesion_multiplier, 0.1, 1e-6);


		double homotypic_multiplier = 4.0;
		force.SetHomotypicPerichondrialSpringConstantMultiplier(
				homotypic_multiplier);

		double alpha = 12.5;
		force.SetAlpha(alpha);
		TS_ASSERT_DELTA(force.GetAlpha(), alpha, 1e-6);

		double repulsion_spring_stiffness = 1.0;
		force.SetRepulsionSpringStiffness(repulsion_spring_stiffness);
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(),
				repulsion_spring_stiffness, 1e-6);

		// Initialise a vector of node forces
		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the x-axis and calculate the force exerted on a neighbour
		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_old_coords =
				cell_population.GetLocationOfCellCentre(my_cell);
		ChastePoint<3> new_coords;
		new_coords.rGetLocation()[0] = my_old_coords[0] + 0.25;
		new_coords.rGetLocation()[1] = my_old_coords[1];
		new_coords.rGetLocation()[2] = my_old_coords[2];
		unsigned node_index = cell_population.GetLocationIndexUsingCell(
				my_cell);
		cell_population.SetNode(node_index, new_coords);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], 0.25,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[2], 0.0,
				1e-4);

		double spring_stiffness_adhesion = force.GetMeinekeSpringStiffness();
		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution =
				force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				spring_stiffness_repulsion * log(0.75), 1e-4); //repulsion force, hence no multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				-0.25 * homotypic_multiplier * spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier, the latter should be one
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 3
		force_contribution = force.CalculateForceBetweenNodes(0, 3,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2],
				0.25 * homotypic_multiplier * baseline_adhesion_multiplier * spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier, the latter should be the baseline_adhesion_multiplier

		// Calculate the force between nodes 0 and 4
		force_contribution = force.CalculateForceBetweenNodes(0, 4,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				(sqrt(2.0) - 1.0) / sqrt(2.0) * homotypic_multiplier
						* ((1.0 - 1.0 / sqrt(2.0)) * (1.0-baseline_adhesion_multiplier) + baseline_adhesion_multiplier)
						* spring_stiffness_adhesion
						* exp(-alpha * (sqrt(2.0) - 1.0)), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2],
				(sqrt(2.0) - 1.0) / sqrt(2.0) * homotypic_multiplier
						* ((1.0 - 1.0/sqrt(2.0)) * 0.9 + 0.1)
						* spring_stiffness_adhesion
						* exp(-alpha * (sqrt(2.0) - 1.0)), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier

		cell_population.Update(false);
		std::vector<std::pair<Node<3>*, Node<3>*> >& r_node_pairs =
				cell_population.rGetNodePairs();
		TS_ASSERT_EQUALS(r_node_pairs.size(), 10);
//
//		// Initialise a vector of node forces
//		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
//			cell_population.GetNode(i)->ClearAppliedForce();
//		}
//		force.AddForceContribution(cell_population);
//
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0],
//				spring_stiffness_repulsion * log(0.75)
//						- 0.25 * homotypic_multiplier
//								* spring_stiffness_adhesion
//								* exp(-alpha * 0.25), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0],
//				0.25 * homotypic_multiplier * spring_stiffness_adhesion
//						* exp(-alpha * 0.25), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0],
//				-spring_stiffness_repulsion * log(0.75), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0,
//				1e-4);

	}

	void TestDirectionalAdhesionForceChondrocyteCase() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, -1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u, false, 1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(3u, false, 0.25, 0.0, 1.25));
		nodes.push_back(new Node<3>(4u, false, 1.25, 0.0, 1.0));
		//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),
				p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		boost::shared_ptr<AbstractCellProperty> p_chondrocyte(
				cell_population.GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>());

		// all cells are chondrocytes
		cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(
				p_chondrocyte);
		cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(
				p_chondrocyte);
		cell_population.GetCellUsingLocationIndex(2)->AddCellProperty(
				p_chondrocyte);
		cell_population.GetCellUsingLocationIndex(3)->AddCellProperty(
				p_chondrocyte);
		cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(
				p_chondrocyte);

		// Create force
		DirectionalAdhesionForce force;

		// Test default values
		double baseline_adhesion_multiplier = force.GetBaselineAdhesionMultiplier();
		TS_ASSERT_DELTA(baseline_adhesion_multiplier, 0.1, 1e-6);


		double homotypic_multiplier = 4.0;
		force.SetHomotypicChondrocyteSpringConstantMultiplier(
				homotypic_multiplier);

		double alpha = 12.5;
		force.SetAlpha(alpha);
		TS_ASSERT_DELTA(force.GetAlpha(), alpha, 1e-6);

		double repulsion_spring_stiffness = 1.0;
		force.SetRepulsionSpringStiffness(repulsion_spring_stiffness);
		TS_ASSERT_DELTA(force.GetRepulsionSpringStiffness(),
				repulsion_spring_stiffness, 1e-6);

		// Initialise a vector of node forces
		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the x-axis and calculate the force exerted on a neighbour
		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_old_coords =
				cell_population.GetLocationOfCellCentre(my_cell);
		ChastePoint<3> new_coords;
		new_coords.rGetLocation()[0] = my_old_coords[0] + 0.25;
		new_coords.rGetLocation()[1] = my_old_coords[1];
		new_coords.rGetLocation()[2] = my_old_coords[2];
		unsigned node_index = cell_population.GetLocationIndexUsingCell(
				my_cell);
		cell_population.SetNode(node_index, new_coords);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], 0.25,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], 0.0,
				1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[2], 0.0,
				1e-4);

		double spring_stiffness_adhesion = force.GetMeinekeSpringStiffness();
		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution =
				force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				spring_stiffness_repulsion * log(0.75), 1e-4); //repulsion force, hence no multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				-0.25 * homotypic_multiplier * baseline_adhesion_multiplier* spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier, the latter should be the baseline adhesion multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 3
		force_contribution = force.CalculateForceBetweenNodes(0, 3,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2],
				0.25 * homotypic_multiplier * spring_stiffness_adhesion
						* exp(-alpha * 0.25), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier, the latter should be one.

		// Calculate the force between nodes 0 and 4
		force_contribution = force.CalculateForceBetweenNodes(0, 4,
				cell_population);
		for (unsigned j = 0; j < 3; j++) {
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0],
				(sqrt(2.0) - 1.0) / sqrt(2.0) * homotypic_multiplier
						* (1.0 / sqrt(2.0) * (1.0-baseline_adhesion_multiplier) + baseline_adhesion_multiplier)
						* spring_stiffness_adhesion
						* exp(-alpha * (sqrt(2.0) - 1.0)), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2],
				(sqrt(2.0) - 1.0) / sqrt(2.0) * homotypic_multiplier
						* (1.0/sqrt(2.0) * (1.0-baseline_adhesion_multiplier) + baseline_adhesion_multiplier)
						* spring_stiffness_adhesion
						* exp(-alpha * (sqrt(2.0) - 1.0)), 1e-4); //adhesion force, hence homotypic multiplier and directional_multiplier

		cell_population.Update(false);
		std::vector<std::pair<Node<3>*, Node<3>*> >& r_node_pairs =
				cell_population.rGetNodePairs();
		TS_ASSERT_EQUALS(r_node_pairs.size(), 10);
//
//		// Initialise a vector of node forces
//		for (unsigned i = 0; i < cell_population.GetNumNodes(); i++) {
//			cell_population.GetNode(i)->ClearAppliedForce();
//		}
//		force.AddForceContribution(cell_population);
//
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0],
//				spring_stiffness_repulsion * log(0.75)
//						- 0.25 * homotypic_multiplier
//								* spring_stiffness_adhesion
//								* exp(-alpha * 0.25), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0],
//				0.25 * homotypic_multiplier * spring_stiffness_adhesion
//						* exp(-alpha * 0.25), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0],
//				-spring_stiffness_repulsion * log(0.75), 1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0,
//				1e-4);
//		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0,
//				1e-4);

	}

	void TestDirectionalAdhesionForceOutputParameters() {
		EXIT_IF_PARALLEL;
		std::string output_directory = "TestForcesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with DirectionalAdhesionForce
		DirectionalAdhesionForce linear_force;
		linear_force.SetCutOffLength(1.5);
		TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "DirectionalAdhesionForce");

		out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("directional_adhesion_results.parameters");
		linear_force.OutputForceParameters(linear_force_parameter_file);
		linear_force_parameter_file->close();

		{
			FileFinder generated_file = output_file_handler.FindFile("directional_adhesion_results.parameters");
			FileFinder reference_file("projects/scaling_cartilage_sheets/test/data/TestDirectionalAdhesionForce/directional_adhesion_results.parameters",
					RelativeTo::ChasteSourceRoot);
			FileComparison comparer(generated_file,reference_file);
			TS_ASSERT(comparer.CompareFiles());
		}

	}

	void TestDirectionalAdhesionForceArchiving() throw (Exception)
	{
		EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DirectionalAdhesionForce.arch";

		{
			DirectionalAdhesionForce force;

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Set member variables
			force.SetMeinekeSpringStiffness(12.34);
			force.SetMeinekeDivisionRestingSpringLength(0.856);
			force.SetMeinekeSpringGrowthDuration(2.593);
			force.SetHomotypicPerichondrialSpringConstantMultiplier(0.051);
			force.SetHomotypicChondrocyteSpringConstantMultiplier(0.091);
			force.SetHeterotypicSpringConstantMultiplier(1.348);
			force.SetBaselineAdhesionMultiplier(0.33);
			force.SetAlpha(12.5);

			// Serialize via pointer to most abstract class possible
			AbstractForce<3>* const p_force = &force;
			output_arch << p_force;
		}

		{
			AbstractForce<3>* p_force;

			// Create an input archive
			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			// Restore from the archive
			input_arch >> p_force;

			// Test member variables
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetHomotypicPerichondrialSpringConstantMultiplier(), 0.051, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetHomotypicChondrocyteSpringConstantMultiplier(), 0.091, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetHeterotypicSpringConstantMultiplier(), 1.348, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalAdhesionForce*>(p_force))->GetAlpha(), 12.5, 1e-6);

			// Tidy up
			delete p_force;
		}
	}

};

#endif /* TESTDIRECTIONALADHESIONFORCE_HPP_ */
