/*
 * NodeBasedCartilageSheet.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: kubuntu1404
 */

#include "NodeBasedCartilageSheet.hpp"

NodeBasedCartilageSheet::NodeBasedCartilageSheet() {
	// TODO Auto-generated constructor stub

}

NodeBasedCartilageSheet::~NodeBasedCartilageSheet() {
	// TODO Auto-generated destructor stub
}

//boost::shared_ptr<NodeBasedCellPopulation<3> > NodeBasedCartilageSheet::Setup() {
//
//	std::vector<Node<3>*> nodes;
//	nodes.push_back(new Node<3>(0u, false, 0.5, 0.0, 0.0));
////	nodes.push_back(new Node<3>(1u, false, -0.5, 0.0, 0.0));
////	nodes.push_back(new Node<3>(2u, false, 0.0, 0.5, 0.0));
////	nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
//
////	NodesOnlyMesh<3> mesh;
//	mesh.ConstructNodesWithoutMesh(nodes, 1.5);
//
////	std::vector<CellPtr> cells;
//	MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
//	CellsGenerator<StochasticDurationCellCycleModel, 3> cells_generator;
//	cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),
//			p_transit_type);
//
//	boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population(new NodeBasedCellPopulation<3>(mesh, cells));
//
//	return cell_population;
//}

boost::shared_ptr<NodeBasedCellPopulation<3> >  NodeBasedCartilageSheet::Setup() {

	bool random_seed = false;
	unsigned n_cells_wide = 3;
	unsigned n_cells_deep = 3;
	unsigned n_cells_high = 3;
	unsigned n_differentiated_cells_width = 0;
	unsigned n_differentiated_cells_depth = 0;
	bool random_birth_times = true;

	/** The next line is needed because this not designed to be run in parallel */
	//EXIT_IF_PARALLEL;
	// TODO change this to a warning and set the problematic input to 1 by default.
	if (n_cells_wide < 1 || n_cells_deep < 1 || n_cells_high < 1) {
		EXCEPTION(
				"The number of cells in x, y or z direction is smaller than 1.");
	}

	std::stringstream ss;

	//unsigned n_cells_per_layer = n_cells_wide*n_cells_deep;
	unsigned n_cells_total = n_cells_wide * n_cells_deep * n_cells_high;

	// Reseed the number generator so that different runs will actually produce different results
	if (random_seed) {
		unsigned seed = time(NULL);
		RandomNumberGenerator::Instance()->Reseed(seed);
		ss << seed;
	}
	std::string filenameaddon_str = ss.str();

	double max_noise = 0.1;

	std::vector<Node<3>*> nodes;
	GenerateRandomHCPNodes(nodes, n_cells_wide, n_cells_deep, n_cells_high,
			max_noise);

	mesh.ConstructNodesWithoutMesh(nodes, 1.5);

	MAKE_PTR(StemCellProliferativeType, p_stem_type);
	MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
	MAKE_PTR(WildTypeCellMutationState, p_state);
	CellsGenerator<CellTissueTypeBasedCellCycleModel, 3> cells_generator;
	cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

	CellTissueTypeBasedCellCycleModel* p_cell_cycle_model =
			new CellTissueTypeBasedCellCycleModel;
	//p_cell_cycle_model->SetDimension(3);

	boost::shared_ptr<NodeBasedCellPopulation<3> > mrCellPopulation(new NodeBasedCellPopulation<3>(mesh, cells));
	//cell_population.SetCellAncestorsToLocationIndices();

	boost::shared_ptr<AbstractCellProperty> p_perichondrial(
			mrCellPopulation->GetCellPropertyRegistry()->Get<
					PerichondrialCellTissueType>());
	boost::shared_ptr<AbstractCellProperty> p_chondrocyte(
			mrCellPopulation->GetCellPropertyRegistry()->Get<
					ChondrocyteCellTissueType>());
	boost::shared_ptr<AbstractCellProperty> p_upwards(
			mrCellPopulation->GetCellPropertyRegistry()->Get<
					UpwardsCellDivisionDirection<3> >());
	boost::shared_ptr<AbstractCellProperty> p_downwards(
			mrCellPopulation->GetCellPropertyRegistry()->Get<
					DownwardsCellDivisionDirection<3> >());

	unsigned lower_index_bottom = n_differentiated_cells_width
			+ n_differentiated_cells_depth * n_cells_wide;
	unsigned upper_index_bottom = n_cells_wide * n_cells_deep
			- n_differentiated_cells_depth * n_cells_wide
			- n_differentiated_cells_width;
	unsigned upper_index_top = n_cells_total
			- (n_differentiated_cells_width
					+ n_differentiated_cells_depth * n_cells_wide);
	unsigned lower_index_top = n_cells_total
			- (n_cells_wide * n_cells_deep
					- n_differentiated_cells_depth * n_cells_wide
					- n_differentiated_cells_width);

	for (AbstractCellPopulation<3>::Iterator cell_iter =
			mrCellPopulation->Begin(); cell_iter != mrCellPopulation->End();
			++cell_iter) {

		unsigned node_index = mrCellPopulation->GetLocationIndexUsingCell(
				*cell_iter);

		// the cell's node_index can be calculated from its coordinates i,j,k (row, column, layer) as
		// node_index = i + j* n_cells_wide + k*n_cells_wide *n_cells_deep,
		// hence reversely they can be obtained via
		unsigned layer_index = node_index / (n_cells_wide * n_cells_deep);//k = node_index/(n_cells_wide *n_cells_deep)  integer division!
		unsigned layer_local_index = node_index % (n_cells_wide * n_cells_deep);// i + j* n_cells_wide = node_index % (n_cells_wide *n_cells_deep)
		unsigned row_index = layer_local_index % n_cells_wide;// i = (i + j* n_cells_wide) % n_cells_wide

		// set the cell tissue type based on the layer index
		if (layer_index == 0 || layer_index == (n_cells_high - 1)) {
			cell_iter->AddCellProperty(p_perichondrial);
		} else {
			cell_iter->AddCellProperty(p_chondrocyte);
		}

		// this gives us the following bounds on the index
		// (if we want stem cells only in the boundary layer and
		// with a padding of a fixed number of differentiated
		// cells in the x and y directions)
		bool stem_cell_bottom = node_index >= lower_index_bottom
				&& node_index < upper_index_bottom
				&& (row_index >= n_differentiated_cells_width)
				&& (row_index < n_cells_wide - n_differentiated_cells_width);
		bool stem_cell_top = node_index >= lower_index_top
				&& node_index < upper_index_top
				&& (row_index >= n_differentiated_cells_width)
				&& (row_index < n_cells_wide - n_differentiated_cells_width);
		if (stem_cell_bottom || stem_cell_top) {
			cell_iter->SetCellProliferativeType(p_stem_type);
			MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (node_index));
			cell_iter->SetAncestor(p_cell_ancestor);

			// set cell division direction
			if (layer_index == 0) {
				cell_iter->AddCellProperty(p_upwards);
			} else {
				cell_iter->AddCellProperty(p_downwards);
			}

			// set random birth times if required
			if (random_birth_times) {
				double birth_time =
						-p_cell_cycle_model->GetAverageStemCellCycleTime()
								* RandomNumberGenerator::Instance()->ranf();
				cell_iter->SetBirthTime(birth_time);
			} else {
				cell_iter->SetBirthTime(-20.0); //Average stem cell cycle time is 24.0 with default values
												//Now we don't have to wait forever for cell divisions to start
			}
		}
	}
	mrCellPopulation->AddCellWriter<CellAncestorWriter>();
	mrCellPopulation->AddCellWriter<CellDivisionDirectionsWriter>();
	mrCellPopulation->AddCellWriter<CellTissueTypesWriter>();

	for (unsigned i = 0; i < nodes.size(); i++) {
		delete nodes[i];
	}

	return mrCellPopulation;
}

/**
 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
 * Arrangement of the nodes will be on a hcp lattice.
 */
void NodeBasedCartilageSheet::GenerateRandomHCPNodes(
		std::vector<Node<3>*> & rNodes, unsigned n_nodes_width,
		unsigned n_nodes_depth, unsigned n_nodes_height, double max_noise) {

	rNodes.clear();
	unsigned n_nodes = n_nodes_width * n_nodes_depth * n_nodes_height;
	rNodes.reserve(n_nodes);

	unsigned id = 0;

	for (unsigned k = 0; k < n_nodes_height; k++) {
		for (unsigned j = 0; j < n_nodes_depth; j++) {
			for (unsigned i = 0; i < n_nodes_width; i++) {
				/*
				 * Note that to pick a random point on the surface of a sphere, it is incorrect
				 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
				 * [0, pi) respectively, since points picked in this way will be 'bunched' near
				 * the poles.
				 */
				double u = RandomNumberGenerator::Instance()->ranf();
				double v = RandomNumberGenerator::Instance()->ranf();

				double noise = max_noise
						* RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = (2 * i + ((j + k) % 2)) * 0.5
						+ noise * cos(random_azimuth_angle)
								* sin(random_zenith_angle);
				double y_coordinate = (sqrt(3) * (j + (k % 2) / 3.0)) * 0.5
						+ noise * sin(random_azimuth_angle)
								* sin(random_zenith_angle);
				double z_coordinate = (2 * sqrt(6) * k / 3.0) * 0.5
						+ noise * cos(random_zenith_angle);
				rNodes.push_back(
						new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
			}
		}
	}
}

