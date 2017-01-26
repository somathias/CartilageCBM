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
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
//#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
//#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
#include "FakePetscSetup.hpp"
#include "OffLatticeSimulation3dDirectedDivision.hpp"
/**
 * Second try to build a model for the cartilage sheet based on a center-based model.
 */

class TestScaling3dNodeBasedCartilageSheets : public AbstractCellBasedTestSuite
{
public:
  
  /**
   * Tests the 3D node-based cartilage sheet simulation
   */
  void Test3dNodeBasedCartilageSheet() throw(Exception)
  {

    bool random_seed = true;
    unsigned n_cells_wide = 10;
    unsigned n_cells_deep = 10;
    unsigned n_cells_high = 1;
    unsigned n_differentiated_cells_width = 0;
    unsigned n_differentiated_cells_depth = 0;
    bool random_birth_times = true;
    
    double spring_stiffness = 1.0;
    double upper_boundary_plane = 4.0;
    double simulation_endtime = 40.0;
    
    std::string output_directory = "3dNodeBasedCartilageSheet/TestGenerateNodesWithAncestors/";
    
    Setup3dNodeBasedCartilageSheet(random_seed, 
				   n_cells_wide,
				   n_cells_deep,
				   n_cells_high,
				   n_differentiated_cells_width,
				   n_differentiated_cells_depth,
				   random_birth_times,
				   spring_stiffness,
				   upper_boundary_plane,
				   output_directory,
				   simulation_endtime );
    
  }

  /**
   * Minimal testing for the generation of the node coordinates
   */
  void Test3dNodeGeneration() throw(Exception)
  {
    unsigned n_nodes_width = 3;
    unsigned n_nodes_depth = 2;
    unsigned n_nodes_height = 1;
    
    std::vector<Node<3>*> nodes;
    GenerateNodes(nodes, n_nodes_width,n_nodes_depth,n_nodes_height);
     
    TS_ASSERT_EQUALS(nodes.size(),n_nodes_width*n_nodes_depth*n_nodes_height );    
    
    for(unsigned i = 0; i< nodes.size(); i++)
    {
      c_vector<double, 3> coordinates = nodes[i]->rGetLocation();
      
      TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.0);
      TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.0);
      TS_ASSERT_EQUALS(coordinates[2], 0.0);
    }
  }
  
private:
  
  /**
   * Sets up the simulation for the specified input parameters
   */
  void Setup3dNodeBasedCartilageSheet(bool random_seed, 
				    unsigned n_cells_wide,
				    unsigned n_cells_deep,
				    unsigned n_cells_high,
				    unsigned n_differentiated_cells_width,
				    unsigned n_differentiated_cells_depth,
				    bool random_birth_times,
				    double spring_stiffness,
				    double upper_boundary_plane,
				    std::string output_directory,
				    double simulation_endtime ) throw(Exception)
  {
    /** The next line is needed because this not designed to be run in parallel */
    EXIT_IF_PARALLEL;
    
    CellBasedEventHandler::Enable();
    std::stringstream ss;
    ss << n_cells_wide << "/";
    ss << n_cells_deep << "/";
    ss << n_cells_high << "/";
    ss << upper_boundary_plane << "/";
    
    //unsigned n_cells_per_layer = n_cells_wide*n_cells_deep;

    
    // Reseed the number generator so that different runs will actually produce different results
    if (random_seed)
    {
      unsigned seed = time(NULL);
      RandomNumberGenerator::Instance()->Reseed(seed); 
      ss << seed;
    }
    std::string filenameaddon_str = ss.str();
    
    std::vector<Node<3>*> nodes;   
    GenerateNodes(nodes, n_cells_wide,n_cells_deep,n_cells_high);
    
    NodesOnlyMesh<3> mesh;
    mesh.ConstructNodesWithoutMesh(nodes, 1.5);
    
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type); 
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); 
    MAKE_PTR(WildTypeCellMutationState, p_state); 
    CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 3> cells_generator;
    cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
     

//     StochasticDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new StochasticDurationGenerationBasedCellCycleModel;
//     p_cell_cycle_model->SetDimension(3);
//     // we could set maxTransitGenerations here.
//     //p_cell_cycle_model->SetMaxTransitGenerations(4);
    
//     
//     //layer of differentiated and stem cells
//     for (unsigned j=0; j<n_cells_deep; j++)
//     {
//       for (unsigned i=0; i<n_cells_wide; i++)
//       {
// 	CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
// 	
// // 	// padding of differentiated cells
// // 	if(i<n_differentiated_cells_width || i>=n_cells_wide-n_differentiated_cells_width || j<n_differentiated_cells_depth || i>=n_cells_deep-n_differentiated_cells_depth)
// // 	{
// // 	  p_cell->SetCellProliferativeType(p_diff_type);
// // 	}
// // 	else
// // 	{
// 	  p_cell->SetCellProliferativeType(p_stem_type);
// 	  MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (i+j*n_cells_wide));
// 	  p_cell->SetAncestor(p_cell_ancestor);
// 	  if(random_birth_times)
// 	  {
// 	    double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
// 	    p_cell->SetBirthTime(birth_time);
// 	  }
// 	  else
// 	  {
// 	    p_cell->SetBirthTime(-20.0); //Average stem cell cycle time is 24.0 with default values 
//                                       //Now we don't have to wait forever for cell divisions to start
// 	  }
// 	  
// // 	}
// 	cells.push_back(p_cell);
//       }
//     }
    
    // n_cells_high-1 more layers of differentiated cells
//     std::vector<CellPtr> cells_extra_layers;
//     cells_generator.GenerateBasicRandom(cells_extra_layers, (n_cells_high-1)*n_cells_per_layer, p_diff_type); 
//     cells.insert(cells.end(),cells_extra_layers.begin(),cells_extra_layers.end());

    NodeBasedCellPopulation<3> cell_population(mesh, cells); 
    cell_population.SetCellAncestorsToLocationIndices();
    cell_population.AddCellWriter<CellAncestorWriter>();
    

    OffLatticeSimulation3dDirectedDivision simulator(cell_population);
    //OffLatticeSimulation<3> simulator(cell_population);
    simulator.SetOutputDirectory(output_directory+filenameaddon_str);
    simulator.SetEndTime(simulation_endtime); //hours
    simulator.SetSamplingTimestepMultiple(12);

    MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
    //p_force->SetCutOffLength(1.5);
    p_force->SetMeinekeSpringStiffness(spring_stiffness);
    simulator.AddForce(p_force);
    
    
    //bottom plane
    c_vector<double,3> point = zero_vector<double>(3);
    c_vector<double,3> normal = zero_vector<double>(3);
    normal(2) = -1.0;
    MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc, (&cell_population, point, normal));
    p_bc->SetUseJiggledNodesOnPlane(true);
    simulator.AddCellPopulationBoundaryCondition(p_bc);
//     
//     //upper plane
//     point(2) = upper_boundary_plane;
//     normal(2) = 1.0;
//     MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_bc_up, (&cell_population, point, normal));
//     p_bc_up->SetUseJiggledNodesOnPlane(true);
//     simulator.AddCellPopulationBoundaryCondition(p_bc_up);

    CellBasedEventHandler::Reset();
    simulator.Solve();
    
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();
    
    for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
  }
  
  /**
   * Generates nodes for a 3D cell sheet a given number of cells wide, deep and high.
   */
  void GenerateNodes(std::vector<Node<3>*> & rNodes, 
		     unsigned n_nodes_width, 
		     unsigned n_nodes_depth,
		     unsigned n_nodes_height) throw(Exception)
  {
    
    rNodes.clear();
    unsigned n_nodes = n_nodes_width*n_nodes_depth*n_nodes_height;
    rNodes.reserve(n_nodes);
    
    unsigned id = 0;
    
    for(unsigned k=0; k<n_nodes_height; k++)
    {
      for(unsigned j=0; j<n_nodes_depth; j++)
      {
	for(unsigned i=0; i<n_nodes_width; i++)
	{
	  rNodes.push_back(new Node<3>(id,  false,  (double) i, (double) j, (double) k));
	  id++;
	}
      }
    }
  }
  
};