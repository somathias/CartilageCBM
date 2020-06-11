/*
 * TestDirectionalRepulsionPWQForce.hpp
 *
 *  Created on: June 10, 2020
 *      Author: Sonja Mathias
 */

#ifndef TESTDIRECTIONALREPULSIONPWQFORCE_HPP_
#define TESTDIRECTIONALREPULSIONPWQFORCE_HPP_

#include <cxxtest/TestSuite.h>
//#include <iostream>

#include "CheckpointArchiveTypes.hpp"

#include "DirectionalRepulsionPWQForce.hpp"
#include "PWQGeneralisedLinearSpringForce.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestDirectionalRepulsionPWQForce: public AbstractCellBasedTestSuite {
public:

	
    void xTestDirectionalRepulsionPWQForcePerpendicular() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, 1.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(2u, false, 0.0, 0.0, 1.0));


		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		// Create force and set adhesion to zero
		DirectionalRepulsionPWQForce force;
        force.SetCutOffLength(1.5);
		force.SetMeinekeSpringStiffness(0.0);

		// Initialise a vector of node forces
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		// Move a node along the z-axis and calculate the force exerted on a neighbour
		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_old_coords = cell_population.GetLocationOfCellCentre(my_cell);
		ChastePoint<3> new_coords;
		new_coords.rGetLocation()[0] = my_old_coords[0];
		new_coords.rGetLocation()[1] = my_old_coords[1];
		new_coords.rGetLocation()[2] = my_old_coords[2]+0.25;
		unsigned node_index = cell_population.GetLocationIndexUsingCell(my_cell);
		cell_population.SetNode(node_index, new_coords);

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetLocation()[2], 0.25, 1e-4);

		double spring_stiffness_adhesion = force.GetMeinekeSpringStiffness();
		TS_ASSERT_DELTA(spring_stiffness_adhesion, 0.0, 1e-4);

		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();
        double cut_off_length = force.GetCutOffLength();

        double ratio = spring_stiffness_adhesion / spring_stiffness_repulsion;
		TS_ASSERT_DELTA(ratio, 0.0, 1e-4);

        double rest_length = 1.0;
        double cut_off_repulsion = rest_length / (1 - sqrt(ratio) * (1-rest_length / cut_off_length));
		TS_ASSERT_DELTA(cut_off_repulsion, 1.0, 1e-4);



		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution = force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

        if (cut_off_repulsion > 0.75)  
        {
            TS_ASSERT_DELTA(force_contribution[2],  
                - spring_stiffness_repulsion*(0.25)*(0.25),
                1e-4);
        }
        else
        {
            TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);
        }  	

		TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

        TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
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

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 
                        - spring_stiffness_repulsion*(0.25)*(0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 
                        spring_stiffness_repulsion*(0.25)*(0.25), 1e-4);

	}

    void TestDirectionalRepulsionPWQForceParallel() {

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.0, 0.0, 0.0));
		nodes.push_back(new Node<3>(1u, false, 0.0, 0.0, 1.0));
		nodes.push_back(new Node<3>(2u, false, 1.0, 0.0, 0.0));


		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		// Create force and set adhesion to zero
		DirectionalRepulsionPWQForce force;
        force.SetCutOffLength(1.5);
		force.SetMeinekeSpringStiffness(0.0);

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
		TS_ASSERT_DELTA(spring_stiffness_adhesion, 0.0, 1e-4);

		double spring_stiffness_repulsion = force.GetRepulsionSpringStiffness();
        double cut_off_length = force.GetCutOffLength();

        double ratio = spring_stiffness_adhesion / spring_stiffness_repulsion;
		TS_ASSERT_DELTA(ratio, 0.0, 1e-4);

        double rest_length = 1.0;
        double cut_off_repulsion = rest_length / (1 - sqrt(ratio) * (1-rest_length / cut_off_length));
		TS_ASSERT_DELTA(cut_off_repulsion, 1.0, 1e-4);


        double baseline_repulsion = force.GetBaselineRepulsionMultiplier();
        TS_ASSERT_EQUALS(baseline_repulsion, 0.1)

		// Calculate the force between nodes 0 and 2
		c_vector<double, 3> force_contribution = force.CalculateForceBetweenNodes(0, 2, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

        if (cut_off_repulsion > 0.75)  
        {
            TS_ASSERT_DELTA(force_contribution[0],  
                - baseline_repulsion*spring_stiffness_repulsion*(0.25)*(0.25),
                1e-4);
        }
        else
        {
            TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
        }  	

		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

        TS_ASSERT_DELTA(force_contribution[0], 0.0, 1e-4);
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

		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 
                        - baseline_repulsion*spring_stiffness_repulsion*(0.25)*(0.25), 
                        1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 
                        baseline_repulsion*spring_stiffness_repulsion*(0.25)*(0.25), 
                        1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0, 1e-4);

	}



	void TestDirectionalRepulsionPWQForceOutputParameters() {
		EXIT_IF_PARALLEL;
		std::string output_directory = "TestForcesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with DirectionalRepulsionPWQForce
		DirectionalRepulsionPWQForce linear_force;
		linear_force.SetCutOffLength(1.5);
		TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "DirectionalRepulsionPWQForce");

		out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("directional_repulsion_results.parameters");
		linear_force.OutputForceParameters(linear_force_parameter_file);
		linear_force_parameter_file->close();

		{
			FileFinder generated_file = output_file_handler.FindFile("directional_repulsion_results.parameters");
			FileFinder reference_file("projects/cartilage/test/data/TestDirectionalRepulsionPWQForce/directional_repulsion_results.parameters",
					RelativeTo::ChasteSourceRoot);
			FileComparison comparer(generated_file,reference_file);
			TS_ASSERT(comparer.CompareFiles());
		}

	}

	void TestDirectionalRepulsionPWQForceArchiving() 
	{
		EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DirectionalRepulsionPWQForce.arch";

		{
			DirectionalRepulsionPWQForce force;

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Set member variables
			force.SetMeinekeSpringStiffness(12.34);
			force.SetMeinekeDivisionRestingSpringLength(0.856);
			force.SetMeinekeSpringGrowthDuration(2.593);
			force.SetBaselineRepulsionMultiplier(0.33);

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
			TS_ASSERT_DELTA((static_cast<DirectionalRepulsionPWQForce*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalRepulsionPWQForce*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalRepulsionPWQForce*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalRepulsionPWQForce*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<DirectionalRepulsionPWQForce*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);

			// Tidy up
			delete p_force;
		}
	}

};

#endif /* TESTDirectionalRepulsionPWQFORCE_HPP_ */
