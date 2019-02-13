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
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
//#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellAncestor.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
#include "FakePetscSetup.hpp"
#include "OffLatticeSimulation2dDirectedDivision.hpp"
/**
 * Second try to build a model for the cartilage sheet based on a center-based model.
 */

class TestScalingNodeBasedCartilageSheets : public AbstractCellBasedTestSuite
{
public:
  
  void TestNodeBasedCartilageSheet() 
  {

    bool random_seed = true;
    unsigned n_layers = 1;
    unsigned n_cells_per_layer = 5;
    unsigned n_differentiated_cells_on_side = 0;
    bool random_birth_times = true;
    
    double spring_stiffness = 1.0;
    double upper_boundary_plane = 4.0;
    double simulation_endtime = 40.0;
    
    std::string output_directory = "NodeBasedCartilageSheet/TestNewSetup/";
    
    SetupNodeBasedCartilageSheet(random_seed, 
				 n_cells_per_layer,
				 n_layers,
				 n_differentiated_cells_on_side,
				 random_birth_times,
				 spring_stiffness,
				 upper_boundary_plane,
				 output_directory,
				 simulation_endtime );
    
  }

private:
  void SetupNodeBasedCartilageSheet(bool random_seed, 
				    unsigned n_cells_per_layer,
				    unsigned n_layers,
				    unsigned n_differentiated_cells_on_side,
				    bool random_birth_times,
				    double spring_stiffness,
				    double upper_boundary_plane,
				    std::string output_directory,
				    double simulation_endtime )
  {
    CellBasedEventHandler::Enable();
    std::stringstream ss;
    ss << n_cells_per_layer << "/";
    ss << n_layers << "/";
    ss << upper_boundary_plane << "/";

    
    // Reseed the number generator so that different runs will actually produce different results
    if (random_seed)
    {
      unsigned seed = time(NULL);
      RandomNumberGenerator::Instance()->Reseed(seed); 
      ss << seed;
    }
    std::string filenameaddon_str = ss.str();
    
    HoneycombMeshGenerator generator(n_cells_per_layer, n_layers); 
    MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
    std::vector<unsigned> location_indices = generator.GetCellLocationIndices(); 
    
    NodesOnlyMesh<2> mesh;
    mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
    
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type); 
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); 
    MAKE_PTR(WildTypeCellMutationState, p_state); 
    CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
    
    //unsigned n_cells = location_indices.size();
    std::vector<CellPtr> cells_current_layer;

    //layer of differentiated and stem cells
    for (unsigned i=0; i<n_cells_per_layer; i++)
    {
      UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel;
      // we could set maxTransitGenerations here.
      //p_cell_cycle_model->SetMaxTransitGenerations(4);
      CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
      if(i<n_differentiated_cells_on_side || i>=n_cells_per_layer-n_differentiated_cells_on_side)
      {
	p_cell->SetCellProliferativeType(p_diff_type);
      }
      else
      {
	p_cell->SetCellProliferativeType(p_stem_type);
	MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (i));
        p_cell->SetAncestor(p_cell_ancestor);
	if(random_birth_times)
	{
	  double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
	  p_cell->SetBirthTime(birth_time);
	}
	else
        {
          p_cell->SetBirthTime(-20.0); //Average stem cell cycle time is 24.0 with default values 
                                      //Now we don't have to wait forever for cell divisions to start
        }


      }
      //p_cell->InitialiseCellCycleModel();
      cells.push_back(p_cell);
    }
    
    // five more layers of differentiated cells
    cells_generator.GenerateBasicRandom(cells_current_layer, (n_layers-1)*n_cells_per_layer, p_diff_type); 
    cells.insert(cells.end(),cells_current_layer.begin(),cells_current_layer.end());

    NodeBasedCellPopulation<2> cell_population(mesh, cells); 
    cell_population.AddCellWriter<CellAncestorWriter>();

    OffLatticeSimulation2dDirectedDivision simulator(cell_population);
    simulator.SetOutputDirectory(output_directory+filenameaddon_str);
    simulator.SetEndTime(simulation_endtime); //hours
    simulator.SetSamplingTimestepMultiple(12);

    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
    //p_force->SetCutOffLength(1.5);
    p_force->SetMeinekeSpringStiffness(spring_stiffness);
    simulator.AddForce(p_force);
    
    
    //bottom plane
    c_vector<double,2> point = zero_vector<double>(2);
    c_vector<double,2> normal = zero_vector<double>(2);
    normal(1) = -1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
    p_bc->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc);
    
    //upper plane
    point(1) = upper_boundary_plane;
    normal(1) = 1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_up, (&cell_population, point, normal));
    p_bc_up->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc_up);

    CellBasedEventHandler::Reset();
    simulator.Solve();
    
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();
  }
  
};