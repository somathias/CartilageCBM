#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FakePetscSetup.hpp"
#include <iostream>

#include "OrientationBasedDivisionRule.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "HorizontalCellDivisionDirection.hpp"
#include "PerturbedUpwardsCellDivisionDirection.hpp"


/**
 * Testing if the orientation based division rule to include directed cell division worked.
 */

class TestOrientationBasedDivisionRule : public AbstractCellBasedTestSuite
{
  public:
    /**
	 * If no CellDivisionDirection is selected, cells should divide randomly.
	 */
    void TestUnDirectedDivision()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        CellPtr my_cell = cell_population.GetCellUsingLocationIndex(0);
        c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new OrientationBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);    

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule = cell_population.GetCentreBasedDivisionRule();
        RandomNumberGenerator::Instance()->Reseed(0);
        std::pair<c_vector<double, 3>, c_vector<double, 3>> positions = p_division_rule->CalculateCellDivisionVector(my_cell, cell_population);

        c_vector<double, 3> your_coords = positions.second;

        c_vector<double, 3> my_new_coords = positions.first;

        RandomNumberGenerator::Instance()->Reseed(0);
        double u = RandomNumberGenerator::Instance()->ranf();
        double v = RandomNumberGenerator::Instance()->ranf();

        double random_azimuth_angle = 2 * M_PI * u;
        double random_zenith_angle = std::acos(2 * v - 1);

        double separation = static_cast<AbstractCentreBasedCellPopulation<3> *>(&(cell_population))->GetMeinekeDivisionSeparation();

        TS_ASSERT_EQUALS(my_coords(0) + 0.5 * separation * cos(random_azimuth_angle) * sin(random_zenith_angle), your_coords(0));
        TS_ASSERT_EQUALS(my_coords(1) + 0.5 * separation * sin(random_azimuth_angle) * sin(random_zenith_angle), your_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) + 0.5 * separation * cos(random_zenith_angle), your_coords(2));

        TS_ASSERT_EQUALS(my_coords(0) - 0.5 * separation * cos(random_azimuth_angle) * sin(random_zenith_angle), my_new_coords(0));
        TS_ASSERT_EQUALS(my_coords(1) - 0.5 * separation * sin(random_azimuth_angle) * sin(random_zenith_angle), my_new_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) - 0.5 * separation * cos(random_zenith_angle), my_new_coords(2));
    }

    void TestUpwardsDirectedDivision()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        MAKE_PTR(UpwardsCellDivisionDirection<3>, p_up_direction);
        CellPtr my_cell = cell_population.GetCellUsingLocationIndex(0);
        my_cell->AddCellProperty(p_up_direction);

        c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

        double separation = static_cast<AbstractCentreBasedCellPopulation<3> *>(&(cell_population))->GetMeinekeDivisionSeparation();

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new OrientationBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);  

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule = cell_population.GetCentreBasedDivisionRule();
        std::pair<c_vector<double, 3>, c_vector<double, 3>> positions = p_division_rule->CalculateCellDivisionVector(my_cell, cell_population);

        c_vector<double, 3> your_coords = positions.second;

        c_vector<double, 3> my_new_coords = positions.first;

        TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
        TS_ASSERT_EQUALS(my_coords(1), your_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) + 0.5 * separation, your_coords(2));

        TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
        TS_ASSERT_EQUALS(my_coords(1), my_new_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) - 0.5 * separation, my_new_coords(2));
    }

    void TestDownwardsDirectedDivision()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        MAKE_PTR(DownwardsCellDivisionDirection<3>, p_down_direction);
        CellPtr my_cell = cell_population.rGetCells().front();
        my_cell->AddCellProperty(p_down_direction);

        c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

        double separation = static_cast<AbstractCentreBasedCellPopulation<3> *>(&(cell_population))->GetMeinekeDivisionSeparation();

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new OrientationBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);  

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule = cell_population.GetCentreBasedDivisionRule();
        std::pair<c_vector<double, 3>, c_vector<double, 3>> positions = p_division_rule->CalculateCellDivisionVector(my_cell, cell_population);

        c_vector<double, 3> your_coords = positions.second;

        c_vector<double, 3> my_new_coords = positions.first;

        TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
        TS_ASSERT_EQUALS(my_coords(1), your_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) - 0.5 * separation, your_coords(2));

        TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
        TS_ASSERT_EQUALS(my_coords(1), my_new_coords(1));
        TS_ASSERT_EQUALS(my_coords(2) + 0.5 * separation, my_new_coords(2));
    }

     void TestHorizontalDirectedDivision()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        MAKE_PTR(HorizontalCellDivisionDirection<3>, p_horizontal_direction);
        CellPtr my_cell = cell_population.rGetCells().front();
        my_cell->AddCellProperty(p_horizontal_direction);

        c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

        double separation = static_cast<AbstractCentreBasedCellPopulation<3> *>(&(cell_population))->GetMeinekeDivisionSeparation();

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new OrientationBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);  

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule = cell_population.GetCentreBasedDivisionRule();
        RandomNumberGenerator::Instance()->Reseed(0);
        std::pair<c_vector<double, 3>, c_vector<double, 3>> positions = p_division_rule->CalculateCellDivisionVector(my_cell, cell_population);

        c_vector<double, 3> your_coords = positions.second;

        c_vector<double, 3> my_new_coords = positions.first;

        RandomNumberGenerator::Instance()->Reseed(0);
        double random_angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

        TS_ASSERT_DELTA(my_coords(0) + 0.5 * separation * cos(random_angle), your_coords(0), 1e-4);
        TS_ASSERT_DELTA(my_coords(1) + 0.5 * separation * sin(random_angle), your_coords(1), 1e-4);
        TS_ASSERT_DELTA(my_coords(2), your_coords(2), 1e-4);

        TS_ASSERT_DELTA(my_coords(0) - 0.5 * separation * cos(random_angle), my_new_coords(0), 1e-4);
        TS_ASSERT_DELTA(my_coords(1) - 0.5 * separation * sin(random_angle), my_new_coords(1), 1e-4);
        TS_ASSERT_DELTA(my_coords(2), my_new_coords(2), 1e-4);
    }
    
    void TestPerturbedUpwardsDirectedDivision()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        std::vector<Node<3> *> nodes;
        nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        //     nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        //     nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        MAKE_PTR(PerturbedUpwardsCellDivisionDirection<3>, p_up_direction);
        p_up_direction->setMaximumZenithAngle(M_PI/8.0);
        CellPtr my_cell = cell_population.GetCellUsingLocationIndex(0);
        my_cell->AddCellProperty(p_up_direction);

        c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);

        double separation = static_cast<AbstractCentreBasedCellPopulation<3> *>(&(cell_population))->GetMeinekeDivisionSeparation();

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new OrientationBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);  

        boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule = cell_population.GetCentreBasedDivisionRule();
        RandomNumberGenerator::Instance()->Reseed(0);

        std::pair<c_vector<double, 3>, c_vector<double, 3>> positions = p_division_rule->CalculateCellDivisionVector(my_cell, cell_population);

        c_vector<double, 3> your_coords = positions.second;

        c_vector<double, 3> my_new_coords = positions.first;
        
        RandomNumberGenerator::Instance()->Reseed(0);
        double u = RandomNumberGenerator::Instance()->ranf();
        double v = RandomNumberGenerator::Instance()->ranf();
        
        double max_zenith_angle = p_up_direction->getMaximumZenithAngle();
        double T = cos(max_zenith_angle);
        double random_azimuth_angle = 2 * M_PI * u;
        double random_zenith_angle = std::acos(T + v * (1 - T));

        TS_ASSERT_DELTA(my_coords(0) + 0.5 * separation * cos(random_azimuth_angle) * sin(random_zenith_angle), your_coords(0), 1e-4);
        TS_ASSERT_DELTA(my_coords(1) + 0.5 * separation * sin(random_azimuth_angle) * sin(random_zenith_angle), your_coords(1), 1e-4);
        TS_ASSERT_DELTA(my_coords(2) + 0.5 * separation * cos(random_zenith_angle), your_coords(2), 1e-4);

        TS_ASSERT_DELTA(my_coords(0) - 0.5 * separation * cos(random_azimuth_angle) * sin(random_zenith_angle), my_new_coords(0), 1e-4);
        TS_ASSERT_DELTA(my_coords(1) - 0.5 * separation * sin(random_azimuth_angle) * sin(random_zenith_angle), my_new_coords(1), 1e-4);
        TS_ASSERT_DELTA(my_coords(2) - 0.5 * separation * cos(random_zenith_angle), my_new_coords(2), 1e-4);
    }

    
};
