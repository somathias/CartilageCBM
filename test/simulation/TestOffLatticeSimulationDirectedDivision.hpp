#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FakePetscSetup.hpp"
#include <iostream>

#include "OffLatticeSimulationDirectedDivision.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"

/**
 * Testing if extending the OffLatticeSimulation class to include directed cell division worked.
 */

class TestOffLatticeSimulationDirectedDivision: public AbstractCellBasedTestSuite {
public:

	/**
	 * If no CellDivisionDirection is selected, cells should divide randomly.
	 */
	void TestUnDirectedDivision() throw (Exception) {
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;



		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
//     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
//     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		CellPtr my_cell = cell_population.rGetCells().front();
		c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

		double separation = static_cast<AbstractCentreBasedCellPopulation<3>*>(&(cell_population))->GetMeinekeDivisionSeparation();

		OffLatticeSimulationDirectedDivision<3> simulator(cell_population);
		RandomNumberGenerator::Instance()->Reseed(0);
		c_vector<double, 3> your_coords = simulator.CalculateCellDivisionVector(my_cell);

		c_vector<double, 3> my_new_coords = cell_population.GetLocationOfCellCentre(my_cell);

		RandomNumberGenerator::Instance()->Reseed(0);
		double u = RandomNumberGenerator::Instance()->ranf();
		double v = RandomNumberGenerator::Instance()->ranf();

		double random_azimuth_angle = 2 * M_PI * u;
		double random_zenith_angle = std::acos(2 * v - 1);

		TS_ASSERT_EQUALS(my_coords(0)+0.5 * separation * cos(random_azimuth_angle)
				* sin(random_zenith_angle), your_coords(0));
		TS_ASSERT_EQUALS(my_coords(1)+ 0.5 * separation * sin(random_azimuth_angle)
				* sin(random_zenith_angle), your_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)+0.5 * separation * cos(random_zenith_angle), your_coords(2));

		TS_ASSERT_EQUALS(my_coords(0)-0.5 * separation * cos(random_azimuth_angle)
				* sin(random_zenith_angle), my_new_coords(0));
		TS_ASSERT_EQUALS(my_coords(1)-0.5 * separation * sin(random_azimuth_angle)
				* sin(random_zenith_angle), my_new_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)-0.5 * separation * cos(random_zenith_angle), my_new_coords(2));

	}

	void TestUpwardsDirectedDivision() throw(Exception)
	{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
		//     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
		//     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
		//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		MAKE_PTR(UpwardsCellDivisionDirection<3>, p_up_direction);
		CellPtr my_cell = cell_population.rGetCells().front();
		my_cell->AddCellProperty(p_up_direction);

		c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

		double separation = static_cast<AbstractCentreBasedCellPopulation<3>*>(&(cell_population))->GetMeinekeDivisionSeparation();

		OffLatticeSimulationDirectedDivision<3> simulator(cell_population);
		c_vector<double, 3> your_coords = simulator.CalculateCellDivisionVector(my_cell);

		c_vector<double, 3> my_new_coords = cell_population.GetLocationOfCellCentre(my_cell);

		TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
		TS_ASSERT_EQUALS(my_coords(1), your_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)+0.5*separation, your_coords(2));

		TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
		TS_ASSERT_EQUALS(my_coords(1), my_new_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)-0.5*separation, my_new_coords(2));

	}

	void TestDownwardsDirectedDivision() throw(Exception)
	{
		/** The next line is needed because we cannot currently run node based simulations in parallel. */
		EXIT_IF_PARALLEL;

		std::vector<Node<3>*> nodes;
		nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
//     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
//     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
//     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);

		MAKE_PTR(DownwardsCellDivisionDirection<3>, p_down_direction);
		CellPtr my_cell = cell_population.rGetCells().front();
		my_cell->AddCellProperty(p_down_direction);

		c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

		double separation = static_cast<AbstractCentreBasedCellPopulation<3>*>(&(cell_population))->GetMeinekeDivisionSeparation();

		OffLatticeSimulationDirectedDivision<3> simulator(cell_population);
		c_vector<double, 3> your_coords = simulator.CalculateCellDivisionVector(my_cell);

		c_vector<double, 3> my_new_coords = cell_population.GetLocationOfCellCentre(my_cell);

		TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
		TS_ASSERT_EQUALS(my_coords(1), your_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)-0.5*separation, your_coords(2));

		TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
		TS_ASSERT_EQUALS(my_coords(1), my_new_coords(1));
		TS_ASSERT_EQUALS(my_coords(2)+0.5*separation, my_new_coords(2));

	}
};
