#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
//#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
//#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "OffLatticeSimulation.hpp"
//#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
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
//     CylindricalHoneycombMeshGenerator generator(20, 1, 2); 
//     Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
    HoneycombMeshGenerator generator(40, 1, 1);    //cells across, up, layers of ghosts
    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    std::vector<unsigned> location_indices = generator.GetCellLocationIndices(); 
    
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type); 
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); 
    MAKE_PTR(WildTypeCellMutationState, p_state); 
    for (unsigned i=0; i<location_indices.size(); i++)
    {
      StochasticDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new  StochasticDurationGenerationBasedCellCycleModel;
      // we could set maxTransitGenerations here.
      CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
      p_cell->SetCellProliferativeType(p_diff_type);
      p_cell->InitialiseCellCycleModel();
      cells.push_back(p_cell);
    }
//     CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
//     cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_stem_type); 
    


    MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); 
    cell_population.SetCellAncestorsToLocationIndices();    
    cell_population.AddPopulationWriter<VoronoiDataWriter>();
    cell_population.AddCellWriter<CellAncestorWriter>();
    cell_population.AddCellWriter<CellAgesWriter>();


    //OffLatticeSimulation<2> simulator(cell_population);
    OffLatticeSimulation2dDirectedDivision simulator(cell_population);
    simulator.SetOutputDirectory("MeshBasedCartilageSheetMaturationManualInitConfig");
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