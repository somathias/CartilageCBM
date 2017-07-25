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
//#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
//#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
//#include "GeneralisedLinearSpringForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellAncestorWriter.hpp"

#include "OffLatticeSimulationDirectedDivision.hpp"
#include "CellTissueTypeBasedCellCycleModel.hpp"
#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"
#include "NodeBasedCartilageSheet.hpp"

#include "PetscSetupAndFinalize.hpp"
/**
 * Second try to build a model for the cartilage sheet based on a center-based model.
 */

class TestScaling3dNodeBasedCartilageSheets: public AbstractCellBasedTestSuite {
public:

	/**
	 * Tests the 3D node-based cartilage sheet simulation
	 */
	void Test3dNodeBasedCartilageSheet() throw (Exception) {

		bool random_seed = false;
		unsigned n_cells_wide = 3;
		unsigned n_cells_deep = 3;
		unsigned n_cells_high = 3;
		unsigned n_differentiated_cells_width = 0;
		unsigned n_differentiated_cells_depth = 0;
		bool random_birth_times = true;

		double spring_stiffness = 1.0;
		double perichondrial_spring_constant_multiplier = 1.0;
		double chondrocyte_spring_constant_multiplier = 1.0;
		double heterotypic_spring_constant_multiplier = 1.0;
		double simulation_endtime = 20.0;

		std::string output_directory =
				"3dNodeBasedCartilageSheet/Test3dHCPCoordinates/";

		RunNodeBasedCartilageSheet(random_seed, n_cells_wide, n_cells_deep,
				n_cells_high, n_differentiated_cells_width,
				n_differentiated_cells_depth, random_birth_times,
				spring_stiffness, perichondrial_spring_constant_multiplier,
				chondrocyte_spring_constant_multiplier,
				heterotypic_spring_constant_multiplier, output_directory,
				simulation_endtime);

	}

	/**
	 * Minimal testing for the generation of the node coordinates
	 */
	void xTest3dNodeGeneration() throw (Exception) {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 2;
		unsigned n_nodes_height = 1;

		std::vector<Node<3>*> nodes;
		GenerateNodes(nodes, n_nodes_width, n_nodes_depth, n_nodes_height);

		TS_ASSERT_EQUALS(nodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		for (unsigned i = 0; i < nodes.size(); i++) {
			c_vector<double, 3> coordinates = nodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.0);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.0);
			TS_ASSERT_EQUALS(coordinates[2], 0.0);
		}
	}

	/**
	 * Minimal testing for the generation of the random node coordinates
	 */
	void xTest3dRandomNodeGeneration() throw (Exception) {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 2;
		unsigned n_nodes_height = 1;
		double max_noise = 0.5;

		std::vector<Node<3>*> nodes;
		GenerateRandomNodes(nodes, n_nodes_width, n_nodes_depth, n_nodes_height,
				max_noise);

		TS_ASSERT_EQUALS(nodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		for (unsigned i = 0; i < nodes.size(); i++) {
			c_vector<double, 3> coordinates = nodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[2], 0.5);
		}
	}

	/**
	 * Minimal testing for the generation of the node coordinates
	 */
	void xTest3dHCPNodeGeneration() throw (Exception) {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 3;
		unsigned n_nodes_height = 3;

		std::vector<Node<3>*> nodes;
		GenerateHCPNodes(nodes, n_nodes_width, n_nodes_depth, n_nodes_height);

		TS_ASSERT_EQUALS(nodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		c_vector<double, 3> coordinates_first = nodes[0]->rGetLocation();
		c_vector<double, 3> coordinates_second = nodes[1]->rGetLocation();
		c_vector<double, 3> coordinates_last =
				nodes[nodes.size() - 1]->rGetLocation();

		//check that first node is in origin
		TS_ASSERT_EQUALS(coordinates_first[0], 0);
		TS_ASSERT_EQUALS(coordinates_first[1], 0);
		TS_ASSERT_EQUALS(coordinates_first[2], 0);

		//check that first two have distance 1 cell diameter
		double distance =
				sqrt(
						(coordinates_first[0] - coordinates_second[0])
								* (coordinates_first[0] - coordinates_second[0])
								+ (coordinates_first[1] - coordinates_second[1])
										* (coordinates_first[1]
												- coordinates_second[1])
								+ (coordinates_first[2] - coordinates_second[2])
										* (coordinates_first[2]
												- coordinates_second[2]));
		TS_ASSERT_DELTA(distance, 1, 1e-4);

		//check the coordinates of the last node
		TS_ASSERT_EQUALS(coordinates_last[0], 2);
		TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-4);
		TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-4);

	}

//	/**
//	 * Minimal testing for the generation of the node coordinates
//	 */
//	void Test3dRandomHCPNodeGeneration() throw (Exception) {
//		unsigned n_nodes_width = 3;
//		unsigned n_nodes_depth = 3;
//		unsigned n_nodes_height = 3;
//
//		std::vector<Node<3>*> nodes;
//		GenerateRandomHCPNodes(nodes, n_nodes_width, n_nodes_depth, n_nodes_height, 0.1);
//
//		TS_ASSERT_EQUALS(nodes.size(),
//				n_nodes_width * n_nodes_depth * n_nodes_height);
//
//		c_vector<double, 3> coordinates_first = nodes[0]->rGetLocation();
//		c_vector<double, 3> coordinates_second = nodes[1]->rGetLocation();
//		c_vector<double, 3> coordinates_last =
//				nodes[nodes.size() - 1]->rGetLocation();
//
//		//check that first node is in origin (+perturbation)
//		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[0], 0.1);
//		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[1], 0.1);
//		TS_ASSERT_LESS_THAN_EQUALS(coordinates_first[2], 0.1);
//
//		//check that first two have distance 1 cell diameter
//		double distance =
//				sqrt(
//						(coordinates_first[0] - coordinates_second[0])
//								* (coordinates_first[0] - coordinates_second[0])
//								+ (coordinates_first[1] - coordinates_second[1])
//										* (coordinates_first[1]
//												- coordinates_second[1])
//								+ (coordinates_first[2] - coordinates_second[2])
//										* (coordinates_first[2]
//												- coordinates_second[2]));
//		TS_ASSERT_LESS_THAN_EQUALS(abs(distance-1), 0.2);
//
//		//check the coordinates of the last node
//		TS_ASSERT_LESS_THAN_EQUALS(coordinates_last[0], 2.1);
//		TS_ASSERT_DELTA(coordinates_last[1], 1.7320, 1e-1);
//		TS_ASSERT_DELTA(coordinates_last[2], 1.6329, 1e-1);
//
//	}

private:

	/**
	 * Sets up the simulation for the specified input parameters
	 */
	void RunNodeBasedCartilageSheet(bool random_seed, unsigned n_cells_wide,
			unsigned n_cells_deep, unsigned n_cells_high,
			unsigned n_differentiated_cells_width,
			unsigned n_differentiated_cells_depth, bool random_birth_times,
			double spring_stiffness,
			double perichondrial_spring_constant_multiplier,
			double chondrocyte_spring_constant_multiplier,
			double heterotypic_spring_constant_multiplier,
			std::string output_directory, double simulation_endtime)
					throw (Exception) {
		CellBasedEventHandler::Enable();

		NodeBasedCartilageSheet* p_cartilage_sheet = new NodeBasedCartilageSheet();

		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population = p_cartilage_sheet->Setup();

		OffLatticeSimulationDirectedDivision<3> simulator(*cell_population);
		//OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory(output_directory);
		simulator.SetEndTime(simulation_endtime);//hours
		simulator.SetSamplingTimestepMultiple(12);

		MAKE_PTR(CellTissueTypeBasedGeneralisedLinearSpringForce<3>, p_force);
		p_force->SetCutOffLength(1.5);
		p_force->SetMeinekeSpringStiffness(spring_stiffness);
		p_force->SetHomotypicPerichondrialSpringConstantMultiplier(perichondrial_spring_constant_multiplier);
		p_force->SetHomotypicChondrocyteSpringConstantMultiplier(chondrocyte_spring_constant_multiplier);
		p_force->SetHeterotypicSpringConstantMultiplier(heterotypic_spring_constant_multiplier);
		simulator.AddForce(p_force);

		CellBasedEventHandler::Reset();
		simulator.Solve();

		CellBasedEventHandler::Headings();
		CellBasedEventHandler::Report();


	}

	/**
	 * Generates nodes for a 3D cell sheet a given number of cells wide, deep and high.
	 * Arrangement of the nodes will be on a Cartesian grid.
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
					rNodes.push_back(new Node<3>(id, false, (double) i, (double) j, (double) k));
					id++;
				}
			}
		}
	}

	/**
	 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
	 * Arrangement of the nodes will be on a Cartesian grid.
	 */
	void GenerateRandomNodes(std::vector<Node<3>*> & rNodes,
			unsigned n_nodes_width,
			unsigned n_nodes_depth,
			unsigned n_nodes_height,
			double max_noise) throw(Exception)
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
					/*
					 * Note that to pick a random point on the surface of a sphere, it is incorrect
					 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
					 * [0, pi) respectively, since points picked in this way will be 'bunched' near
					 * the poles.
					 */
					double u = RandomNumberGenerator::Instance()->ranf();
					double v = RandomNumberGenerator::Instance()->ranf();

					double noise = max_noise * RandomNumberGenerator::Instance()->ranf();

					double random_azimuth_angle = 2 * M_PI * u;
					double random_zenith_angle = std::acos(2 * v - 1);

					double x_coordinate = i + noise * cos(random_azimuth_angle) * sin(random_zenith_angle);
					double y_coordinate = j + noise * sin(random_azimuth_angle) * sin(random_zenith_angle);
					double z_coordinate = k + noise * cos(random_zenith_angle);
					rNodes.push_back(new Node<3>(id, false, x_coordinate, y_coordinate, z_coordinate));
					id++;
				}
			}
		}
	}

	/**
	 * Generates nodes for a 3D cell sheet a given number of cells wide, deep and high.
	 * Arrangement of the nodes will be on a hexagonal close packed (hcp) lattice.
	 */
	void GenerateHCPNodes(std::vector<Node<3>*> & rNodes,
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
					double x_coordinate = (2*i + ((j+k) % 2))*0.5; //cell radius = 0.5 cell diameter (reference length)
					double y_coordinate = (sqrt(3)*(j+(k % 2)/3.0))*0.5;//cell radius = 0.5 cell diameter (reference length)
					double z_coordinate = (2*sqrt(6)*k/3.0)*0.5;//cell radius = 0.5 cell diameter (reference length)
					rNodes.push_back(new Node<3>(id, false, x_coordinate, y_coordinate, z_coordinate));
					id++;
				}
			}
		}
	}


};
