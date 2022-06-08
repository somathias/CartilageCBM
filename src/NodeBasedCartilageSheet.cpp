/*
 * NodeBasedCartilageSheet.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: kubuntu1404
 */

#include "NodeBasedCartilageSheet.hpp"

NodeBasedCartilageSheet::NodeBasedCartilageSheet() : mNumberOfNodesPerXDimension(3), mNumberOfNodesPerYDimension(3),
													 mNumberOfNodesPerZDimension(3), 
													 mNumberOfPerichondrialLayersAbove(0),
													 mNumberOfPerichondrialLayersBelow(1),
													 mMaxCoordinatePerturbation(0), 
													 mStemCellG1Duration(60.0),
													 mTransitCellG1Duration(30.0),
													 mSPhaseDuration(10.0),
													 mSeed(0), mPatchSizeLimit(6), 
													 mSynchronizeCellCycles(false),
													 mDivisionDirections(true),
													 mNodesGenerated(false),
													 mCellPopulationSetup(false)
{
													
													 
}

NodeBasedCartilageSheet::~NodeBasedCartilageSheet()
{
	// TODO Auto-generated destructor stub
}

void NodeBasedCartilageSheet::Setup()
{

	// mesh generation
	if (!mNodesGenerated)
	{
		GenerateNodesOnCartesianGrid(1.0);
	}
	else if (mNodes.size() != mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension * mNumberOfNodesPerZDimension)
	{
		EXCEPTION(
			"Sheet dimensions and number of generated coordinates do not match.");
	}
	mMesh.ConstructNodesWithoutMesh(mNodes, 1.5);

	//nodes themselves can be deleted (the mesh copies them)
	for (unsigned i = 0; i < mNodes.size(); i++)
	{
		delete mNodes[i];
	}
	mNodesGenerated = false;

	//generate the cells
	mCells.clear(); //this code is copied from the CellsGenerator class - not sure if needed
	mCells.reserve(mMesh.GetNumNodes());  //this code is copied from the CellsGenerator class - not sure if needed
	MAKE_PTR(WildTypeCellMutationState, p_state);
	MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
	for (unsigned i=0; i<mMesh.GetNumNodes(); i++)
    {
		CellTissueTypeBasedCellCycleModel* p_model = new CellTissueTypeBasedCellCycleModel();
		p_model->SetStemCellG1Duration(mStemCellG1Duration);
		p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
		p_model->SetSDuration(mSPhaseDuration);
		p_model->SetMDuration(1e-12);
		p_model->SetG2Duration(1e-12);
		p_model->SetMaxTransitGenerations(2);
		p_model->SetPatchSizeLimit(mPatchSizeLimit);

		CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);
		// setting of the birth time will be done later when initialising the random stem cell configuration
		mCells.push_back(p_cell);
	}

	//CellsGenerator<CellTissueTypeBasedCellCycleModel, 3> cells_generator;
	//cells_generator.GenerateBasicRandom(mCells, mMesh.GetNumNodes(), p_diff_type);

	//make copy of mCells because calling reset will clear the vector passed to the cell population
	std::vector<CellPtr> tmp_mCells = mCells;

	//generate the cell population
	mpCellPopulation.reset(new NodeBasedCellPopulation<3>(mMesh, tmp_mCells));


	//set the division rule
	boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule_to_set(new OrientationBasedDivisionRule<3, 3>());
	mpCellPopulation->SetCentreBasedDivisionRule(p_division_rule_to_set);

	mCellPopulationSetup = true;
}

bool NodeBasedCartilageSheet::isCellPopulationSetup() const
{
	return mCellPopulationSetup;
}

boost::shared_ptr<NodeBasedCellPopulation<3>> NodeBasedCartilageSheet::GetCellPopulation()
{
	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}
	return mpCellPopulation;
}

void NodeBasedCartilageSheet::InitialiseTissueLayersAndCellDivisionDirections()
{

	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}

	boost::shared_ptr<AbstractCellProperty> p_perichondrial(
		mpCellPopulation->GetCellPropertyRegistry()->Get<PerichondrialCellTissueType>());
	boost::shared_ptr<AbstractCellProperty> p_chondrocyte(
		mpCellPopulation->GetCellPropertyRegistry()->Get<ChondrocyteCellTissueType>());
	boost::shared_ptr<AbstractCellProperty> p_upwards(
		mpCellPopulation->GetCellPropertyRegistry()->Get<UpwardsCellDivisionDirection<3>>());
	boost::shared_ptr<AbstractCellProperty> p_downwards(
		mpCellPopulation->GetCellPropertyRegistry()->Get<DownwardsCellDivisionDirection<3>>());
	boost::shared_ptr<AbstractCellProperty> p_horizontal(
		mpCellPopulation->GetCellPropertyRegistry()->Get<HorizontalCellDivisionDirection<3>>());
	boost::shared_ptr<AbstractCellProperty> p_upper_layer(
		mpCellPopulation->GetCellPropertyRegistry()->Get<UpperPerichondrialLayer>());
	boost::shared_ptr<AbstractCellProperty> p_lower_layer(
		mpCellPopulation->GetCellPropertyRegistry()->Get<LowerPerichondrialLayer>());

	for (AbstractCellPopulation<3>::Iterator cell_iter =
			 mpCellPopulation->Begin();
		 cell_iter != mpCellPopulation->End();
		 ++cell_iter)
	{
		unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(
			*cell_iter);

		// the cell's node_index can be calculated from its coordinates i,j,k (row, column, layer) as
		// node_index = i + j* n_cells_wide + k*n_cells_wide *n_cells_deep,
		// hence reversely they can be obtained via
		unsigned layer_index = node_index / (mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension); //k = node_index/(n_cells_wide *n_cells_deep)  integer division!
																										 //		unsigned layer_local_index = node_index
																										 //				% (mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension);// i + j* n_cells_wide = node_index % (n_cells_wide *n_cells_deep)
																										 //		unsigned row_index = layer_local_index % mNumberOfNodesPerXDimension;// i = (i + j* n_cells_wide) % n_cells_wide

		// set the cell tissue type and cell division direction based on the layer index
		if (layer_index <= mNumberOfPerichondrialLayersBelow - 1)
		{
			cell_iter->AddCellProperty(p_perichondrial);
			if(mDivisionDirections){
				cell_iter->AddCellProperty(p_horizontal);
				cell_iter->AddCellProperty(p_lower_layer);
			}
			
		}
		else if (layer_index >= mNumberOfNodesPerZDimension - mNumberOfPerichondrialLayersAbove)
		{
			cell_iter->AddCellProperty(p_perichondrial);
			if (mDivisionDirections) {
				cell_iter->AddCellProperty(p_horizontal);
				cell_iter->AddCellProperty(p_upper_layer);
			}
		}
		else
		{
			cell_iter->AddCellProperty(p_chondrocyte);
		}
	}

	mpCellPopulation->AddCellWriter<CellDivisionDirectionsWriter>();
	mpCellPopulation->AddCellWriter<CellTissueTypesWriter>();
}

void NodeBasedCartilageSheet::InitialiseBulkStemCellConfiguration(
	unsigned numberOfCellsWide, unsigned numberOfCellsDeep)
{

	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}
	unsigned n_cells_total = mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension * mNumberOfNodesPerZDimension;
	unsigned n_differentiated_cells_width = (mNumberOfNodesPerXDimension - numberOfCellsWide) / 2;
	unsigned n_differentiated_cells_depth = (mNumberOfNodesPerYDimension - numberOfCellsDeep) / 2;

	unsigned lower_index_bottom = n_differentiated_cells_width + n_differentiated_cells_depth * mNumberOfNodesPerXDimension;
	unsigned upper_index_bottom = mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension - n_differentiated_cells_depth * mNumberOfNodesPerXDimension - n_differentiated_cells_width;
	unsigned upper_index_top =
		n_cells_total - (n_differentiated_cells_width + n_differentiated_cells_depth * mNumberOfNodesPerXDimension);
	unsigned lower_index_top = n_cells_total - (mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension - n_differentiated_cells_depth * mNumberOfNodesPerXDimension - n_differentiated_cells_width);

	MAKE_PTR(StemCellProliferativeType, p_stem_type);

	for (AbstractCellPopulation<3>::Iterator cell_iter =
			 mpCellPopulation->Begin();
		 cell_iter != mpCellPopulation->End();
		 ++cell_iter)
	{
		unsigned node_index = mpCellPopulation->GetLocationIndexUsingCell(
			*cell_iter);

		// the cell's node_index can be calculated from its coordinates i,j,k (row, column, layer) as
		// node_index = i + j* n_cells_wide + k*n_cells_wide *n_cells_deep,
		// hence reversely they can be obtained via
		//		unsigned layer_index = node_index
		//				/ (mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension);//k = node_index/(n_cells_wide *n_cells_deep)  integer division!
		unsigned layer_local_index = node_index % (mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension); // i + j* n_cells_wide = node_index % (n_cells_wide *n_cells_deep)
		unsigned row_index = layer_local_index % mNumberOfNodesPerXDimension;								   // i = (i + j* n_cells_wide) % n_cells_wide

		// this gives us the following bounds on the index
		// (if we want stem cells only in the boundary layer and
		// with a padding of a fixed number of differentiated
		// cells in the x and y directions)
		bool stem_cell_bottom = node_index >= lower_index_bottom && node_index < upper_index_bottom && (row_index >= n_differentiated_cells_width) && (row_index < mNumberOfNodesPerXDimension - n_differentiated_cells_width);
		bool stem_cell_top = node_index >= lower_index_top && node_index < upper_index_top && (row_index >= n_differentiated_cells_width) && (row_index < mNumberOfNodesPerXDimension - n_differentiated_cells_width);

		if (stem_cell_bottom || stem_cell_top)
		{

			cell_iter->SetCellProliferativeType(p_stem_type);
			MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (node_index));
			cell_iter->SetAncestor(p_cell_ancestor);

			// set random birth times if required
			if (!mSynchronizeCellCycles)
			{
				CellTissueTypeBasedCellCycleModel *p_cell_cycle_model =
					new CellTissueTypeBasedCellCycleModel;
				double birth_time =
					-p_cell_cycle_model->GetAverageStemCellCycleTime() * RandomNumberGenerator::Instance()->ranf();
				cell_iter->SetBirthTime(birth_time);
			}
			else
			{
				cell_iter->SetBirthTime(-56.0); //Average stem cell cycle time is 60.0 with current values
												//Now we don't have to wait forever for cell divisions to start
			}
		}
	}
	mpCellPopulation->AddCellWriter<CellAncestorWriter>();
}


void NodeBasedCartilageSheet::InitialiseRandomStemCellConfiguration(
	unsigned numberOfStemCells)
{

	//check if the population is set up
	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}

	// sanity check input parameter
	unsigned n_cells_total = mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension * mNumberOfNodesPerZDimension;
	unsigned n_cells_total_per_layer = mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension;
	if (numberOfStemCells > n_cells_total)
	{
		EXCEPTION(
			"Specified number of stem cells larger than total number of cells.");
	}
	unsigned numberOfStemCellsPerLayer;
	if (mNumberOfNodesPerZDimension == 1)
	{
		numberOfStemCellsPerLayer = numberOfStemCells;
	}
	else
	{
		unsigned numberOfPerichondrialLayers = mNumberOfPerichondrialLayersAbove + mNumberOfPerichondrialLayersBelow;
		numberOfStemCellsPerLayer = floor(numberOfStemCells / static_cast<double>(numberOfPerichondrialLayers));
	}

	if (mNumberOfNodesPerXDimension > 2 && mNumberOfNodesPerYDimension > 2 && numberOfStemCellsPerLayer > (mNumberOfNodesPerXDimension-2)*(mNumberOfNodesPerYDimension-2)){
		EXCEPTION(
			"Specified number of stem cells per layer does not allow for leaving a boundary of non-activated stem cells to prevent boundary effects.");
	}

	MAKE_PTR(StemCellProliferativeType, p_stem_type);

	// generate stem cell indices for all layers
	// This evenly distributes the specified total number of stem cells among all layers of perichondrial cells
	unsigned i = 0;
	while (i < numberOfStemCells)
	{
		//choose a row
		unsigned row;
		if(mNumberOfNodesPerXDimension > 2){
			// if at least 3 rows available leave outer rows untouched
			row = RandomNumberGenerator::Instance()->randMod(mNumberOfNodesPerXDimension-2) + 1;
		}
		else{
			row = RandomNumberGenerator::Instance()->randMod(mNumberOfNodesPerXDimension);
		}
		//choose a column		
		unsigned column;
		if(mNumberOfNodesPerYDimension > 2){
			// if at least 3 columns available leave outer columns untouched
			column = RandomNumberGenerator::Instance()->randMod(mNumberOfNodesPerYDimension-2) + 1;
		}
		else{
			column = RandomNumberGenerator::Instance()->randMod(mNumberOfNodesPerYDimension);
		}

		//calculate node index
		// if in lower layer the offset is zero, else the offset is n_cells_total_per_layer*(mNumberOfNodesPerZDimension -1);
		unsigned offset = (i < numberOfStemCellsPerLayer) ? 0 : n_cells_total_per_layer * (mNumberOfNodesPerZDimension - 1);
		unsigned node_index = row + column * mNumberOfNodesPerXDimension + offset;

		//get cell belonging to node index
		CellPtr cell = mpCellPopulation->GetCellUsingLocationIndex(node_index);

		// set proliferative type to stem cell if differentiated in order to not choose the same cell twice
		if (!cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
		{

			cell->SetCellProliferativeType(p_stem_type);

			MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (node_index));
			cell->SetAncestor(p_cell_ancestor);

			// // set random birth times if required
			// if (!mSynchronizeCellCycles)
			// {
			// 	CellTissueTypeBasedCellCycleModel *p_cell_cycle_model =
			// 		new CellTissueTypeBasedCellCycleModel;
			// 	double birth_time =
			// 		-p_cell_cycle_model->GetAverageStemCellCycleTime() * RandomNumberGenerator::Instance()->ranf();
			// 	cell->SetBirthTime(birth_time);
			// }
			// else
			// {
			// 	cell->SetBirthTime(-56.0); //Average stem cell cycle time is 60.0 with current values
			// 							   //Now we don't have to wait forever for cell divisions to start
			// }
			cell->SetBirthTime(0.0);
			//increase counter
			i++;
		}
	}
	mpCellPopulation->AddCellWriter<CellAncestorWriter>();
}

void NodeBasedCartilageSheet::InitialiseMissingColumnExperiment()
{
	//check if the population is set up
	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}

	// currently hard-code setup for sheet of 5x5xZ cells
	if (mNumberOfNodesPerXDimension != 5)
	{
		EXCEPTION(
			"This experiment setup is currently hard-coded for 5 cells in the x-direction.");
	}
	if (mNumberOfNodesPerYDimension != 5)
	{
		EXCEPTION(
			"This experiment setup is currently hard-coded for 5 cells in the y-direction.");
	}

	//get index of column in the center - currently hard-coded as index 12
	unsigned center_index = 12;
	for (unsigned layer = 0; layer < mNumberOfNodesPerZDimension; layer++)
	{
		unsigned index = center_index + layer*mNumberOfNodesPerXDimension*mNumberOfNodesPerYDimension;
		mMesh.DeleteNode(index);

		//CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(index);
		//mpCellPopulation->RemoveCellUsingLocationIndex(index, p_cell); //is this necessary since I regenerate the cell population anyways.
		mCells.erase(mCells.begin()+index); 
	}

	//generate the cell population again, because I can't delete cells from mpCellPopulation->mCells
	mpCellPopulation.reset(new NodeBasedCellPopulation<3>(mMesh, mCells));

	//set the division rule again
	boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule_to_set(new OrientationBasedDivisionRule<3, 3>());
	mpCellPopulation->SetCentreBasedDivisionRule(p_division_rule_to_set);

	mCellPopulationSetup = true;

	// make cell with location index 6 divide into space provided by deleting the column
	boost::shared_ptr<AbstractCellProperty> p_perichondrial(
		mpCellPopulation->GetCellPropertyRegistry()->Get<PerichondrialCellTissueType>());
	boost::shared_ptr<AbstractCellProperty> p_lower_layer(
		mpCellPopulation->GetCellPropertyRegistry()->Get<LowerPerichondrialLayer>());
	boost::shared_ptr<AbstractCellProperty> p_fixed(
		mpCellPopulation->GetCellPropertyRegistry()->Get<FixedCellDivisionDirection<3>>());

	unsigned index = 6;
	CellPtr p_cell = mpCellPopulation->GetCellUsingLocationIndex(index);

	p_cell->AddCellProperty(p_perichondrial);
	p_cell->AddCellProperty(p_lower_layer);
	p_cell->AddCellProperty(p_fixed);

	MAKE_PTR(StemCellProliferativeType, p_stem_type);
	p_cell->SetCellProliferativeType(p_stem_type);

	MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (index));
	p_cell->SetAncestor(p_cell_ancestor);

	// set random birth times if required
	if (!mSynchronizeCellCycles)
	{
		CellTissueTypeBasedCellCycleModel *p_cell_cycle_model =
			new CellTissueTypeBasedCellCycleModel;
		double birth_time =
			-p_cell_cycle_model->GetAverageStemCellCycleTime() * RandomNumberGenerator::Instance()->ranf();
		p_cell->SetBirthTime(birth_time);
	}
	else
	{
		p_cell->SetBirthTime(-56.0); //Average stem cell cycle time is 60.0 with current values
									//Now we don't have to wait forever for cell divisions to start
	}

	mpCellPopulation->AddCellWriter<CellDivisionDirectionsWriter>();
	mpCellPopulation->AddCellWriter<CellTissueTypesWriter>();
	mpCellPopulation->AddCellWriter<CellAncestorWriter>();			

}


void NodeBasedCartilageSheet::SetCartilageSheetDimensions(
	unsigned numberOfCellsWide, unsigned numberOfCellsDeep,
	unsigned numberOfCellsHigh)
{
	mNumberOfNodesPerXDimension = numberOfCellsWide;
	mNumberOfNodesPerYDimension = numberOfCellsDeep;
	mNumberOfNodesPerZDimension = numberOfCellsHigh;
}

void NodeBasedCartilageSheet::SetPatchSizeLimit(unsigned patchSizeLimit){
	assert(patchSizeLimit >= 1);
	mPatchSizeLimit = patchSizeLimit;
}

void NodeBasedCartilageSheet::SetPhaseDurations(double stem_cell_g1_duration, double transit_cell_g1_duration, double s_phase_duration){
	mStemCellG1Duration = stem_cell_g1_duration;
	mTransitCellG1Duration = transit_cell_g1_duration;
	mSPhaseDuration = s_phase_duration;
}

/**
 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
 * Arrangement of the nodes will be on a Cartesian grid.
 */
void NodeBasedCartilageSheet::GenerateNodesOnCartesianGrid(double scaling)
{

	mNodes.clear();
	unsigned n_nodes_width = mNumberOfNodesPerXDimension;
	unsigned n_nodes_depth = mNumberOfNodesPerYDimension;
	unsigned n_nodes_height = mNumberOfNodesPerZDimension;
	unsigned n_nodes = n_nodes_width * n_nodes_depth * n_nodes_height;
	mNodes.reserve(n_nodes);

	unsigned id = 0;

	for (unsigned k = 0; k < n_nodes_height; k++)
	{
		for (unsigned j = 0; j < n_nodes_depth; j++)
		{
			for (unsigned i = 0; i < n_nodes_width; i++)
			{
				/*
				 * Note that to pick a random point on the surface of a sphere, it is incorrect
				 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
				 * [0, pi) respectively, since points picked in this way will be 'bunched' near
				 * the poles.
				 */
				double u = RandomNumberGenerator::Instance()->ranf();
				double v = RandomNumberGenerator::Instance()->ranf();

				double noise = mMaxCoordinatePerturbation * RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = scaling*(i + noise * cos(random_azimuth_angle) * sin(random_zenith_angle));
				double y_coordinate = scaling*(j + noise * sin(random_azimuth_angle) * sin(random_zenith_angle));
				double z_coordinate = k + noise * cos(random_zenith_angle);
				mNodes.push_back(
					new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
			}
		}
	}
	mNodesGenerated = true;
}

double NodeBasedCartilageSheet::getMaxCoordinatePerturbation() const
{
	return mMaxCoordinatePerturbation;
}

void NodeBasedCartilageSheet::setMaxCoordinatePerturbation(
	double maxCoordinatePerturbation)
{
	mMaxCoordinatePerturbation = maxCoordinatePerturbation;
}

unsigned NodeBasedCartilageSheet::getNumberOfNodesPerXDimension() const
{
	return mNumberOfNodesPerXDimension;
}

unsigned NodeBasedCartilageSheet::getNumberOfNodesPerYDimension() const
{
	return mNumberOfNodesPerYDimension;
}

unsigned NodeBasedCartilageSheet::getNumberOfNodesPerZDimension() const
{
	return mNumberOfNodesPerZDimension;
}

unsigned NodeBasedCartilageSheet::getNumberOfPerichondrialLayersAbove() const
{
	return mNumberOfPerichondrialLayersAbove;
}

unsigned NodeBasedCartilageSheet::getNumberOfPerichondrialLayersBelow() const
{
	return mNumberOfPerichondrialLayersBelow;
}

void NodeBasedCartilageSheet::setNumberOfPerichondrialLayersAbove(unsigned numberOfPerichondrialLayersAbove)
{
	mNumberOfPerichondrialLayersAbove = numberOfPerichondrialLayersAbove;
}

void NodeBasedCartilageSheet::setNumberOfPerichondrialLayersBelow(unsigned numberOfPerichondrialLayersBelow)
{
	mNumberOfPerichondrialLayersBelow = numberOfPerichondrialLayersBelow;
}

unsigned NodeBasedCartilageSheet::getSeed() const
{
	return mSeed;
}

void NodeBasedCartilageSheet::setSynchronizeCellCycles(
	bool synchronizeCellCycles)
{
	mSynchronizeCellCycles = synchronizeCellCycles;
}

void NodeBasedCartilageSheet::setDivisionDirections(
	bool divisionDirections)
{
	mDivisionDirections = divisionDirections;
}

void NodeBasedCartilageSheet::UseRandomSeed()
{

	// Reseed the number generator so that different runs will actually produce different results
	mSeed = time(NULL);
	RandomNumberGenerator::Instance()->Reseed(mSeed);
}

/**
 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
 * Arrangement of the nodes will be on a hcp lattice.
 */
void NodeBasedCartilageSheet::GenerateNodesOnHCPGrid(double scaling)
{

	mNodes.clear();
	unsigned n_nodes_width = mNumberOfNodesPerXDimension;
	unsigned n_nodes_depth = mNumberOfNodesPerYDimension;
	unsigned n_nodes_height = mNumberOfNodesPerZDimension;
	unsigned n_nodes = n_nodes_width * n_nodes_depth * n_nodes_height;
	mNodes.reserve(n_nodes);

	unsigned id = 0;

	for (unsigned k = 0; k < n_nodes_height; k++)
	{
		for (unsigned j = 0; j < n_nodes_depth; j++)
		{
			for (unsigned i = 0; i < n_nodes_width; i++)
			{
				/*
				 * Note that to pick a random point on the surface of a sphere, it is incorrect
				 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
				 * [0, pi) respectively, since points picked in this way will be 'bunched' near
				 * the poles.
				 */
				double u = RandomNumberGenerator::Instance()->ranf();
				double v = RandomNumberGenerator::Instance()->ranf();

				double noise = mMaxCoordinatePerturbation * RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = scaling*((2 * i + ((j + k) % 2)) * 0.5 + noise * cos(random_azimuth_angle) * sin(random_zenith_angle));
				double y_coordinate = scaling*((sqrt(3) * (j + (k % 2) / 3.0)) * 0.5 + noise * sin(random_azimuth_angle) * sin(random_zenith_angle));
				double z_coordinate = (2 * sqrt(6) * k / 3.0) * 0.5 + noise * cos(random_zenith_angle);
				mNodes.push_back(
					new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
			}
		}
	}
	mNodesGenerated = true;
}

void NodeBasedCartilageSheet::GenerateNodesOnStackedHexagonalGrid(double scaling)
{

	mNodes.clear();
	unsigned n_nodes_width = mNumberOfNodesPerXDimension;
	unsigned n_nodes_depth = mNumberOfNodesPerYDimension;
	unsigned n_nodes_height = mNumberOfNodesPerZDimension;
	unsigned n_nodes = n_nodes_width * n_nodes_depth * n_nodes_height;
	mNodes.reserve(n_nodes);

	unsigned id = 0;

	for (unsigned k = 0; k < n_nodes_height; k++)
	{
		for (unsigned j = 0; j < n_nodes_depth; j++)
		{
			for (unsigned i = 0; i < n_nodes_width; i++)
			{
				/*
				 * Note that to pick a random point on the surface of a sphere, it is incorrect
				 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
				 * [0, pi) respectively, since points picked in this way will be 'bunched' near
				 * the poles.
				 */
				double u = RandomNumberGenerator::Instance()->ranf();
				double v = RandomNumberGenerator::Instance()->ranf();

				double noise = mMaxCoordinatePerturbation * RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = scaling*((2 * i + (j % 2)) * 0.5 + noise * cos(random_azimuth_angle) * sin(random_zenith_angle));
				double y_coordinate = scaling*((sqrt(3) * j) * 0.5 + noise * sin(random_azimuth_angle) * sin(random_zenith_angle));
				double z_coordinate = k + noise * cos(random_zenith_angle);
				mNodes.push_back(
					new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
			}
		}
	}
	mNodesGenerated = true;
}

