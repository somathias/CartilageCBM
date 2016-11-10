#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "FakePetscSetup.hpp"
#include "OffLatticeSimulation2dDirectedDivision.hpp"
/**
 * First try to build a model for the cartilage sheet based on a center-based model using a Voronoi tesselation.
 */

class TestScalingMeshBasedCartilageSheets : public AbstractCellBasedTestSuite
{
public:
  void TestMeshBasedCartilageSheet() throw(Exception)
  {
    CylindricalHoneycombMeshGenerator generator(20, 1, 2); 
    Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
//     HoneycombMeshGenerator generator(20, 1, 4);    //cells across, up, layers of ghosts
//     MutableMesh<2,2>* p_mesh = generator.GetMesh();
    std::vector<unsigned> location_indices = generator.GetCellLocationIndices(); 
    
    std::vector<CellPtr> cells;
    MAKE_PTR(TransitCellProliferativeType, p_transit_type); 
    CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
    cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type); 

    MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); 
    cell_population.SetCellAncestorsToLocationIndices();    
    cell_population.AddPopulationWriter<VoronoiDataWriter>();
    cell_population.AddCellWriter<CellAncestorWriter>();


    //OffLatticeSimulation<2> simulator(cell_population);
    OffLatticeSimulation2dDirectedDivision simulator(cell_population);
    simulator.SetOutputDirectory("MeshBasedCartilageSheetPeriodicBC");
    simulator.SetEndTime(50.0); // what unit is this??? Seems to be hours
    simulator.SetSamplingTimestepMultiple(12);

    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
    p_force->SetCutOffLength(1.5);
    simulator.AddForce(p_force);
    
    c_vector<double,2> point = zero_vector<double>(2);
    c_vector<double,2> normal = zero_vector<double>(2);
    normal(1) = -1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
    simulator.AddCellPopulationBoundaryCondition(p_bc);

    simulator.Solve();
  }
};