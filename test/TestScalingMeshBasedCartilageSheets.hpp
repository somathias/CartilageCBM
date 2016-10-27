#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

/**
 * First try to build a model for the cartilage sheet based on a center-based model using a Voronoi tesselation.
 */

class TestScalingMeshBasedCartilageSheets : public AbstractCellBasedTestSuite
{
public:
  void TestMeshBasedCartilageSheet() throw(Exception)
  {
    HoneycombMeshGenerator generator(6, 4, 4);    // Parameters are: cells across, cells up, number of ghost cell layers
    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    std::vector<unsigned> location_indices = generator.GetCellLocationIndices(); //necessary when using ghost nodes
    
    std::vector<CellPtr> cells;
    MAKE_PTR(TransitCellProliferativeType, p_transit_type); // the type influences the average length of the G1 phase, ie. the time between cell divisions
    CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
    //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_stem_type); 
    cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type); //necessary when using ghost nodes
    
    //MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
    MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //necessary when using ghost nodes
    
    cell_population.AddPopulationWriter<VoronoiDataWriter>();

    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetOutputDirectory("MeshBasedCartilageSheetSolidBottomBoundary");
    simulator.SetEndTime(25.0); // what unit is this???
    simulator.SetSamplingTimestepMultiple(12);

    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
    //p_force->SetCutOffLength(1.5);
    simulator.AddForce(p_force);
    
    c_vector<double,2> point = zero_vector<double>(2);
    c_vector<double,2> normal = zero_vector<double>(2);
    normal(1) = -1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
    simulator.AddCellPopulationBoundaryCondition(p_bc);

    simulator.Solve();

    //TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
    //TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
  }
};