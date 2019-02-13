#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "FakePetscSetup.hpp"
#include <iostream>

#include "OffLatticeSimulation2dDirectedDivision.hpp"

/**
 * Testing if extending the OffLatticeSimulation class to include directed cell division worked.
 */

class TestOffLatticeSimulation2dDirectedDivision : public AbstractCellBasedTestSuite
{
public:
  void TestUpwardsDirectedDivision() 
  {
    
    HoneycombMeshGenerator generator(1, 1);
    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    
    std::vector<CellPtr> cells;
    CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

    MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
    
    CellPtr my_cell = cell_population.rGetCells().front();
    c_vector<double, 2> my_coords = cell_population.GetLocationOfCellCentre(my_cell);
    
    double separation = static_cast<AbstractCentreBasedCellPopulation<2>*>(&(cell_population))->GetMeinekeDivisionSeparation();
    
    OffLatticeSimulation2dDirectedDivision simulator(cell_population);
    c_vector<double, 2> your_coords = simulator.CalculateCellDivisionVector(my_cell);
    
    c_vector<double, 2> my_new_coords = cell_population.GetLocationOfCellCentre(my_cell);
    
    TS_ASSERT_EQUALS(my_coords(0), your_coords(0));
    TS_ASSERT_EQUALS(my_coords(1)+0.5*separation, your_coords(1));
    
    TS_ASSERT_EQUALS(my_coords(0), my_new_coords(0));
    TS_ASSERT_EQUALS(my_coords(1)-0.5*separation, my_new_coords(1));
  
  }
};