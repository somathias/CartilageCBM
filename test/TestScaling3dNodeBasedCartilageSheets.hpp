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

#include "FakePetscSetup.hpp"
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
		unsigned n_cells_wide = 2;
		unsigned n_cells_deep = 4;
		unsigned n_cells_high = 1;
		unsigned n_differentiated_cells_width = 1;
		unsigned n_differentiated_cells_depth = 2;
		bool random_birth_times = true;

		double spring_stiffness = 1.0;
		double perichondrial_spring_constant_multiplier = 1.0;
		double chondrocyte_spring_constant_multiplier = 1.0;
		double heterotypic_spring_constant_multiplier = 1.0;
		double simulation_endtime = 40.0;

		std::string output_directory =
				"3dNodeBasedCartilageSheet/Test3dSheetRandomNodeCoordinates/noise-0-1/";

		Setup3dNodeBasedCartilageSheet(random_seed, n_cells_wide, n_cells_deep,
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
	void Test3dNodeGeneration() throw (Exception) {
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
	void Test3dRandomNodeGeneration() throw (Exception) {
		unsigned n_nodes_width = 3;
		unsigned n_nodes_depth = 2;
		unsigned n_nodes_height = 1;
		double max_noise = 0.5;

		std::vector<Node<3>*> nodes;
		GenerateRandomNodes(nodes, n_nodes_width, n_nodes_depth, n_nodes_height, max_noise);

		TS_ASSERT_EQUALS(nodes.size(),
				n_nodes_width * n_nodes_depth * n_nodes_height);

		for (unsigned i = 0; i < nodes.size(); i++) {
			c_vector<double, 3> coordinates = nodes[i]->rGetLocation();

			TS_ASSERT_LESS_THAN_EQUALS(coordinates[0], 2.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[1], 1.5);
			TS_ASSERT_LESS_THAN_EQUALS(coordinates[2], 0.5);
		}
	}

private:

	/**
	 * Sets up the simulation for the specified input parameters
	 */
	void Setup3dNodeBasedCartilageSheet(bool random_seed, unsigned n_cells_wide,
			unsigned n_cells_deep, unsigned n_cells_high,
			unsigned n_differentiated_cells_width,
			unsigned n_differentiated_cells_depth, bool random_birth_times,
			double spring_stiffness,
			double perichondrial_spring_constant_multiplier,
			double chondrocyte_spring_constant_multiplier,
			double heterotypic_spring_constant_multiplier,
			std::string output_directory, double simulation_endtime)
					throw (Exception) {
		/** The next line is needed because this not designed to be run in parallel */
		EXIT_IF_PARALLEL;

		// TODO change this to a warning and set the problematic input to 1 by default.
		if (n_cells_wide <1 || n_cells_deep < 1 || n_cells_high <1 )
		{
			EXCEPTION("The number of cells in x, y or z direction is smaller than 1.");
		}

		CellBasedEventHandler::Enable();
		std::stringstream ss;
		ss << n_cells_wide << "/";
		ss << n_cells_deep << "/";
		ss << n_cells_high << "/";
		ss << spring_stiffness << "/";
		ss << perichondrial_spring_constant_multiplier << "/";
		ss << chondrocyte_spring_constant_multiplier << "/";
		ss << heterotypic_spring_constant_multiplier << "/";

		//unsigned n_cells_per_layer = n_cells_wide*n_cells_deep;
		unsigned n_cells_total = n_cells_wide*n_cells_deep*n_cells_high;

		// Reseed the number generator so that different runs will actually produce different results
		if (random_seed)
		{
			unsigned seed = time(NULL);
			RandomNumberGenerator::Instance()->Reseed(seed);
			ss << seed;
		}
		std::string filenameaddon_str = ss.str();

		double max_noise = 0.1;

		std::vector<Node<3>*> nodes;
		GenerateRandomNodes(nodes, n_cells_wide,n_cells_deep,n_cells_high, max_noise);

		NodesOnlyMesh<3> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(WildTypeCellMutationState, p_state);
		CellsGenerator<CellTissueTypeBasedCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

		CellTissueTypeBasedCellCycleModel* p_cell_cycle_model = new CellTissueTypeBasedCellCycleModel;
		//p_cell_cycle_model->SetDimension(3);

		NodeBasedCellPopulation<3> cell_population(mesh, cells);
		//cell_population.SetCellAncestorsToLocationIndices();

		boost::shared_ptr<AbstractCellProperty> p_perichondrial(cell_population.GetCellPropertyRegistry()->Get<PerichondrialCellTissueType>());
		boost::shared_ptr<AbstractCellProperty> p_chondrocyte(cell_population.GetCellPropertyRegistry()->Get<ChondrocyteCellTissueType>());
		boost::shared_ptr<AbstractCellProperty> p_upwards(cell_population.GetCellPropertyRegistry()->Get<UpwardsCellDivisionDirection<3> >());
		boost::shared_ptr<AbstractCellProperty> p_downwards(cell_population.GetCellPropertyRegistry()->Get<DownwardsCellDivisionDirection<3> >());

		unsigned lower_index_bottom = n_differentiated_cells_width + n_differentiated_cells_depth*n_cells_wide;
		unsigned upper_index_bottom = n_cells_wide*n_cells_deep - n_differentiated_cells_depth*n_cells_wide - n_differentiated_cells_width;
		unsigned upper_index_top = n_cells_total - (n_differentiated_cells_width + n_differentiated_cells_depth*n_cells_wide);
		unsigned lower_index_top = n_cells_total - (n_cells_wide*n_cells_deep - n_differentiated_cells_depth*n_cells_wide - n_differentiated_cells_width);

		for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{

			unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

			// the cell's node_index can be calculated from its coordinates i,j,k (row, column, layer) as
			// node_index = i + j* n_cells_wide + k*n_cells_wide *n_cells_deep,
			// hence reversely they can be obtained via
			unsigned layer_index = node_index/(n_cells_wide *n_cells_deep);//k = node_index/(n_cells_wide *n_cells_deep)  integer division!
			unsigned layer_local_index = node_index % (n_cells_wide*n_cells_deep);// i + j* n_cells_wide = node_index % (n_cells_wide *n_cells_deep)
			unsigned row_index = layer_local_index % n_cells_wide;// i = (i + j* n_cells_wide) % n_cells_wide

			// set the cell tissue type based on the layer index
			if(layer_index == 0 || layer_index == (n_cells_high-1)) {
				cell_iter->AddCellProperty(p_perichondrial);
			} else {
				cell_iter->AddCellProperty(p_chondrocyte);
			}

			// this gives us the following bounds on the index
			// (if we want stem cells only in the boundary layer and
			// with a padding of a fixed number of differentiated
			// cells in the x and y directions)
			bool stem_cell_bottom = node_index >= lower_index_bottom && node_index < upper_index_bottom && (row_index >= n_differentiated_cells_width) && (row_index < n_cells_wide - n_differentiated_cells_width);
			bool stem_cell_top = node_index >= lower_index_top && node_index < upper_index_top && (row_index >= n_differentiated_cells_width) && (row_index < n_cells_wide - n_differentiated_cells_width);
			if (stem_cell_bottom || stem_cell_top)
			{
				cell_iter->SetCellProliferativeType(p_stem_type);
				MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (node_index));
				cell_iter->SetAncestor(p_cell_ancestor);

				// set cell division direction
				if (layer_index == 0) {
					cell_iter->AddCellProperty(p_upwards);
				}
				else {
					cell_iter->AddCellProperty(p_downwards);
				}

				// set random birth times if required
				if(random_birth_times)
				{
					double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
					cell_iter->SetBirthTime(birth_time);
				}
				else
				{
					cell_iter->SetBirthTime(-20.0); //Average stem cell cycle time is 24.0 with default values
													//Now we don't have to wait forever for cell divisions to start
				}
			}
		}
		cell_population.AddCellWriter<CellAncestorWriter>();
		cell_population.AddCellWriter<CellDivisionDirectionsWriter>();
		cell_population.AddCellWriter<CellTissueTypesWriter>();

		OffLatticeSimulationDirectedDivision<3> simulator(cell_population);
		//OffLatticeSimulation<3> simulator(cell_population);
		simulator.SetOutputDirectory(output_directory+filenameaddon_str);
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
					rNodes.push_back(new Node<3>(id, false, (double) i, (double) j, (double) k));
					id++;
				}
			}
		}
	}

	/**
	 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
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
					 * the poles. See #2230.
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

};
