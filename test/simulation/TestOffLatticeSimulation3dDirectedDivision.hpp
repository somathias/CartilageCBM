#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "FakePetscSetup.hpp"
#include <iostream>

#include "OffLatticeSimulation3dDirectedDivision.hpp"

/**
 * Testing if extending the OffLatticeSimulation class to include directed cell division worked.
 */

class TestOffLatticeSimulation3dDirectedDivision : public AbstractCellBasedTestSuite
{
public:
  void TestUpwardsDirectedDivision()
  {
    /** The next line is needed because we cannot currently run node based simulations in parallel. */
    EXIT_IF_PARALLEL;
    
    std::vector<Node<3>*> nodes;
    nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
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
    
    CellPtr my_cell = cell_population.rGetCells().front();
    c_vector<double, 3> my_coords = cell_population.GetLocationOfCellCentre(my_cell);
    
    double separation = static_cast<AbstractCentreBasedCellPopulation<3>*>(&(cell_population))->GetMeinekeDivisionSeparation();
    
    OffLatticeSimulation3dDirectedDivision simulator(cell_population);
    c_vector<double, 3> your_coords = simulator.CalculateCellDivisionVector(my_cell);
    
    c_vector<double, 3> my_new_coords = cell_population.GetLocationOfCellCentre(my_cell);
    
    TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
    TS_ASSERT_EQUALS(my_coords(1), your_coords(1));
    TS_ASSERT_EQUALS(my_coords(2)+0.5*separation, your_coords(2));
    
    TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
    TS_ASSERT_EQUALS(my_coords(1), my_new_coords(1));
    TS_ASSERT_EQUALS(my_coords(2)-0.5*separation, my_new_coords(2));
  
  }
};