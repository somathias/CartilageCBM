#include <ctime>
#include <string>
#include <sstream>
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
//#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "OffLatticeSimulation.hpp"
//#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellAncestor.hpp"
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
  void TestMeshBasedCartilageSheet() 
  {
    CellBasedEventHandler::Enable();
    
    // Reseed the number generator so that different runs will actually produce different results
    unsigned seed = time(NULL);
    std::stringstream ss;
    ss << seed;
    std::string seed_str = ss.str();
    RandomNumberGenerator::Instance()->Reseed(seed); 
//     CylindricalHoneycombMeshGenerator generator(20, 1, 2); 
//     Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
    unsigned n_layers = 6;
    
    HoneycombMeshGenerator generator(15, n_layers, 1);    //cells across, up, layers of ghosts
    MutableMesh<2,2>* p_mesh = generator.GetMesh();
    std::vector<unsigned> location_indices = generator.GetCellLocationIndices(); 
    
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type); 
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); 
    MAKE_PTR(WildTypeCellMutationState, p_state); 
    CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
    
    unsigned n_cells = location_indices.size();
    unsigned n_cells_per_layer = n_cells/n_layers;
  
    //no bottom layers of differentiated cells
    std::vector<CellPtr> cells_current_layer;
//     cells_generator.GenerateBasicRandom(cells_current_layer, 2*n_cells_per_layer, p_diff_type); 
//     cells.insert(cells.end(),cells_current_layer.begin(),cells_current_layer.end());
    
    //layer of differentiated and stem cells
    for (unsigned i=0; i<n_cells_per_layer; i++)
    {
      UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new  UniformG1GenerationalCellCycleModel;
      //p_cell_cycle_model->SetMaxTransitGenerations(4);
      // we could set maxTransitGenerations here.
      CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
      if(i<7 || i>=n_cells_per_layer-7)
      {
	p_cell->SetCellProliferativeType(p_diff_type);
      }
      else
      {
	p_cell->SetCellProliferativeType(p_stem_type);
	double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
	p_cell->SetBirthTime(birth_time);
	MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (i));
        p_cell->SetAncestor(p_cell_ancestor);
      }
      //p_cell->InitialiseCellCycleModel();
      cells.push_back(p_cell);
    }
    
    // five more layers of differentiated cells
    cells_generator.GenerateBasicRandom(cells_current_layer, 5*n_cells_per_layer, p_diff_type); 
    cells.insert(cells.end(),cells_current_layer.begin(),cells_current_layer.end());
    
    //check if we initialised the correct number of cells
    TS_ASSERT_EQUALS(cells.size(), n_cells);


    MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); 
    //MeshBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices); 
    //cell_population.SetCellAncestorsToLocationIndices();    
    cell_population.AddPopulationWriter<VoronoiDataWriter>();
    cell_population.AddCellWriter<CellAncestorWriter>();
    cell_population.AddCellWriter<CellAgesWriter>();


    //OffLatticeSimulation<2> simulator(cell_population);
    OffLatticeSimulation2dDirectedDivision simulator(cell_population);
    simulator.SetOutputDirectory("MeshBasedCartilageSheetSingleStemCell/"+seed_str);
    simulator.SetEndTime(100.0); // what unit is this??? Seems to be hours
    simulator.SetSamplingTimestepMultiple(12);

    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
    //p_force->SetCutOffLength(1.5);
    simulator.AddForce(p_force);
    
    
    //bottom plane
    c_vector<double,2> point = zero_vector<double>(2);
    c_vector<double,2> normal = zero_vector<double>(2);
    normal(1) = -1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
    p_bc->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc);
    
    //upper plane
    point(1) = 4.0;
    normal(1) = 1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_up, (&cell_population, point, normal));
    p_bc_up->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc_up);

    simulator.Solve();
    
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();
  }
};