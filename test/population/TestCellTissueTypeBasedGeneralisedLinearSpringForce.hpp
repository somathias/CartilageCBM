#ifndef TESTCELLTISSUETYPEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_
#define TESTCELLTISSUETYPEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
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

class TestCellTissueTypeBasedGeneralisedLinearSpringForce: public AbstractCellBasedTestSuite {
public:

	void TestCellTissueTypeBasedGeneralisedLinearSpringForceWithMeshBasedCellPopulation()
			throw (Exception) {
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
		CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasic(cells, location_indices.size(), location_indices);

		// Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		// Create force
		CellTissueTypeBasedGeneralisedLinearSpringForce<2> force;

		// Test set/get method
		TS_ASSERT_DELTA(force.GetHomotypicPerichondrialSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHomotypicChondrocyteSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetAlpha(), 5.0, 1e-6);

		force.SetHomotypicPerichondrialSpringConstantMultiplier(2.0);
		force.SetHomotypicChondrocyteSpringConstantMultiplier(3.0);
		force.SetHeterotypicSpringConstantMultiplier(4.0);
//		force.SetAlpha(12.5);

		TS_ASSERT_DELTA(force.GetHomotypicPerichondrialSpringConstantMultiplier(), 2.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHomotypicChondrocyteSpringConstantMultiplier(), 3.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 4.0, 1e-6);
//		TS_ASSERT_DELTA(force.GetAlpha(), 12.5, 1e-6);

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

		// Test the case where node 59 and its neighbours do not have a cell tissue type
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

		// Next, test the case where node 59 is a perichondrial but its neighbours are chondrocytes...
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		boost::shared_ptr<AbstractCellProperty> p_chondrocyte(cell_population.GetCellPropertyRegistry()->Get<ChondrocyteCellTissueType>());
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			cell_iter->AddCellProperty(p_chondrocyte);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<PerichondrialCellTissueType>(), false);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<ChondrocyteCellTissueType>(), true);
		}

		boost::shared_ptr<AbstractCellProperty> p_perichondrial(cell_population.GetCellPropertyRegistry()->Get<PerichondrialCellTissueType>());
		cell_population.GetCellUsingLocationIndex(59)->RemoveCellProperty<ChondrocyteCellTissueType>();
		cell_population.GetCellUsingLocationIndex(59)->AddCellProperty(p_perichondrial);

		TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(59)->HasCellProperty<PerichondrialCellTissueType>(), true);
		TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(59)->HasCellProperty<ChondrocyteCellTissueType>(), false);

		force.AddForceContribution(cell_population);

		// ...for which the force magnitude should be increased by 4, our chosen multiplier for heterotypic interactions under attraction
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 4.0*0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+4.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

		// Next, test the case where node 59 is a chondrocyte but its neighbours are perichondrial...
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			cell_iter->RemoveCellProperty<ChondrocyteCellTissueType>(); //remove cell tissue type from test before
			cell_iter->AddCellProperty(p_perichondrial);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<PerichondrialCellTissueType>(), true);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<ChondrocyteCellTissueType>(), false);
		}

		cell_population.GetCellUsingLocationIndex(59)->RemoveCellProperty<PerichondrialCellTissueType>();
		cell_population.GetCellUsingLocationIndex(59)->AddCellProperty(p_chondrocyte);

		TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(59)->HasCellProperty<PerichondrialCellTissueType>(), false);
		TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(59)->HasCellProperty<ChondrocyteCellTissueType>(), true);

		force.AddForceContribution(cell_population);

		// ...for which the force magnitude should also be increased by 4, our chosen multiplier for heterotypic interactions under attraction
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 4.0*0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+4.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

		// test the case where node 59 and its neighbours are all perichondrial cells...
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		cell_population.GetCellUsingLocationIndex(59)->RemoveCellProperty<ChondrocyteCellTissueType>(); //delete its type from the last test

		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			cell_iter->AddCellProperty(p_perichondrial);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<PerichondrialCellTissueType>(), true);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<ChondrocyteCellTissueType>(), false);
		}

		force.AddForceContribution(cell_population);

		// ...for which the force magnitude should be increased by 2, our chosen multiplier for homotypic perichondrial interactions, again only for attractive interactions
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 2.0*0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+2.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);

		// Finally, test the case where node 59 and its neighbours are all chondrocytes cells...
		for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
		{
			cell_population.GetNode(i)->ClearAppliedForce();
		}

		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			cell_iter->RemoveCellProperty<PerichondrialCellTissueType>(); //remove cell tissue type from test before
			cell_iter->AddCellProperty(p_chondrocyte);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<PerichondrialCellTissueType>(), false);
			TS_ASSERT_EQUALS(cell_iter->HasCellProperty<ChondrocyteCellTissueType>(), true);
		}

		force.AddForceContribution(cell_population);

		// ...for which the force magnitude should be increased by 3, our chosen multiplier for homotypic chondrocyte interactions, again only for attractive interactions
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[0], 3.0*0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(58)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[0], -0.5*spring_stiffness+3.0*(4.0/sqrt(7.0)-2.5)*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(59)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[0], 0.5*spring_stiffness, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(60)->rGetAppliedForce()[1], 0.0, 1e-4);
	}

	void TestCellTissueTypeBasedGeneralisedLinearSpringForceWithNodeBasedCellPopulation() {

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
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		// Create force
		CellTissueTypeBasedGeneralisedLinearSpringForce<3> force;

		// Test set/get method
		TS_ASSERT_DELTA(force.GetHomotypicPerichondrialSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHomotypicChondrocyteSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetHeterotypicSpringConstantMultiplier(), 1.0, 1e-6);
		TS_ASSERT_DELTA(force.GetAlpha(), 5.0, 1e-6);

		double alpha = 12.5;
		force.SetAlpha(alpha);
		TS_ASSERT_DELTA(force.GetAlpha(), alpha, 1e-6);

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
		std::cout << node_index << std::endl;
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

		TS_ASSERT_DELTA(force_contribution[0], spring_stiffness*log(0.75), 1e-4);
		TS_ASSERT_DELTA(force_contribution[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(force_contribution[2], 0.0, 1e-4);

		// Calculate the force between nodes 0 and 1
		force_contribution = force.CalculateForceBetweenNodes(0, 1, cell_population);
		for (unsigned j=0; j<3; j++)
		{
			assert(!std::isnan(force_contribution[j]));
		}

		TS_ASSERT_DELTA(force_contribution[0], -0.25*spring_stiffness*exp(-alpha*0.25), 1e-4);
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



		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], spring_stiffness*log(0.75) - 0.25*spring_stiffness*exp(-alpha*0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], 0.25*spring_stiffness*exp(-alpha*0.25), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], -spring_stiffness*log(0.75), 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
		TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 0.0, 1e-4);

	}
	void TestCellTissueTypeBasedGeneralisedLinearSpringForceOutputParameters()
	{
		EXIT_IF_PARALLEL;
		std::string output_directory = "TestForcesOutputParameters";
		OutputFileHandler output_file_handler(output_directory, false);

		// Test with CellTissueTypeBasedGeneralisedLinearSpringForce
		CellTissueTypeBasedGeneralisedLinearSpringForce<2> linear_force;
		linear_force.SetCutOffLength(1.5);
		TS_ASSERT_EQUALS(linear_force.GetIdentifier(), "CellTissueTypeBasedGeneralisedLinearSpringForce-2-2");

		out_stream linear_force_parameter_file = output_file_handler.OpenOutputFile("tissue_type_based_results.parameters");
		linear_force.OutputForceParameters(linear_force_parameter_file);
		linear_force_parameter_file->close();

		{
			FileFinder generated_file = output_file_handler.FindFile("tissue_type_based_results.parameters");
			FileFinder reference_file("projects/scaling_cartilage_sheets/test/data/TestCellTissueTypeBasedGeneralisedLinearSpringForce/tissue_type_based_results.parameters",
					RelativeTo::ChasteSourceRoot);
			FileComparison comparer(generated_file,reference_file);
			TS_ASSERT(comparer.CompareFiles());
		}

	}

	void TestCellTissueTypeBasedGeneralisedLinearSpringForceArchiving() throw (Exception)
	{
		EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellTissueTypeBasedGeneralisedLinearSpringForce.arch";

		{
			CellTissueTypeBasedGeneralisedLinearSpringForce<2> force;

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			// Set member variables
			force.SetMeinekeSpringStiffness(12.34);
			force.SetMeinekeDivisionRestingSpringLength(0.856);
			force.SetMeinekeSpringGrowthDuration(2.593);
			force.SetHomotypicPerichondrialSpringConstantMultiplier(0.051);
			force.SetHomotypicChondrocyteSpringConstantMultiplier(0.091);
			force.SetHeterotypicSpringConstantMultiplier(1.348);
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
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringStiffness(), 12.34, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeDivisionRestingSpringLength(), 0.856, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetMeinekeSpringGrowthDuration(), 2.593, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetHomotypicPerichondrialSpringConstantMultiplier(), 0.051, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetHomotypicChondrocyteSpringConstantMultiplier(), 0.091, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetHeterotypicSpringConstantMultiplier(), 1.348, 1e-6);
			TS_ASSERT_DELTA((static_cast<CellTissueTypeBasedGeneralisedLinearSpringForce<2>*>(p_force))->GetAlpha(), 12.5, 1e-6);

			// Tidy up
			delete p_force;
		}
	}

};

#endif /*TESTCELLTISSUETYPEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_*/
