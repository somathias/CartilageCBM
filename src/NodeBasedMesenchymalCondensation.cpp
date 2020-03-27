/*
 * NodeBasedMesenchymalCondensation.cpp
 *
 *  Created on: Mar 12, 2020
 *      Author: Sonja Mathias
 */

#include "NodeBasedMesenchymalCondensation.hpp"

NodeBasedMesenchymalCondensation::NodeBasedMesenchymalCondensation() : mNumberOfNodesPerXDimension(3), mNumberOfNodesPerYDimension(3),
													 mMaxCoordinatePerturbation(0), mDistanceBetweeenBoundaries(7.0),
													 mSeed(0), mSynchronizeCellCycles(false),
													 mDivisionDirections(true),
													 mNodesGenerated(false),
													 mCellPopulationSetup(false)
{
}

NodeBasedMesenchymalCondensation::~NodeBasedMesenchymalCondensation()
{
	// TODO Auto-generated destructor stub
}

void NodeBasedMesenchymalCondensation::Setup()
{

	// mesh generation
	if (!mNodesGenerated)
	{
		GenerateNodesOnCartesianGrid();
	}
	else if (mNodes.size() != mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension)
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
		ChondrocytesOnlyCellCycleModel* p_model = new ChondrocytesOnlyCellCycleModel();
		p_model->SetTransitCellG1Duration(10.0);
		p_model->SetSDuration(3.0);
		p_model->SetMDuration(1e-12);
		p_model->SetG2Duration(1e-12);
		p_model->SetMaxTransitGenerations(2);
		CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->SetCellProliferativeType(p_diff_type);

		// setting of the birth time will be done later when initialising the random transit cell configuration
		mCells.push_back(p_cell);
	}

	//generate the cell population
	mpCellPopulation.reset(new NodeBasedCellPopulation<3>(mMesh, mCells));

	//set the division rule
	boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>> p_division_rule_to_set(new OrientationBasedDivisionRule<3, 3>());
	mpCellPopulation->SetCentreBasedDivisionRule(p_division_rule_to_set);

    mpCellPopulation->AddCellWriter<CellDivisionDirectionsWriter>();

	mCellPopulationSetup = true;

}

bool NodeBasedMesenchymalCondensation::isCellPopulationSetup() const
{
	return mCellPopulationSetup;
}

boost::shared_ptr<NodeBasedCellPopulation<3>> NodeBasedMesenchymalCondensation::GetCellPopulation()
{
	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}
	return mpCellPopulation;
}


void NodeBasedMesenchymalCondensation::InitialiseRandomConfiguration(
	unsigned numberOfCells)
{

	//check if the population is set up
	if (!mCellPopulationSetup)
	{
		EXCEPTION("The cell population has not been set up yet.");
	}

	// sanity check input parameter
	unsigned n_cells_total = mNumberOfNodesPerXDimension * mNumberOfNodesPerYDimension;
	if (numberOfCells > n_cells_total)
	{
		EXCEPTION(
			"Specified number of activated cells larger than total number of cells.");
	}

	MAKE_PTR(TransitCellProliferativeType, p_transit_type);
	boost::shared_ptr<AbstractCellProperty> p_upwards(
		mpCellPopulation->GetCellPropertyRegistry()->Get<UpwardsCellDivisionDirection<3>>());


	// generate transit cell indices for all layers
	// This evenly distributes the specified total number of transit cells among all layers of perichondrial cells
	unsigned i = 0;
	while (i < numberOfCells)
	{
		//choose a row
		unsigned row = RandomNumberGenerator::Instance()->randMod(
			mNumberOfNodesPerXDimension);
		//choose a column
		unsigned column = RandomNumberGenerator::Instance()->randMod(
			mNumberOfNodesPerYDimension);
		//calculate node index
		unsigned node_index = row + column * mNumberOfNodesPerXDimension;

		//get cell belonging to node index
		CellPtr cell = mpCellPopulation->GetCellUsingLocationIndex(node_index);

		// set proliferative type to transit cell if differentiated in order to not choose the same cell twice
		if (!cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
		{

			cell->SetCellProliferativeType(p_transit_type);

			MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (node_index));
			cell->SetAncestor(p_cell_ancestor);

			if(mDivisionDirections){
			cell->AddCellProperty(p_upwards);
		}

			// set random birth times if required
			if (!mSynchronizeCellCycles)
			{
				ChondrocytesOnlyCellCycleModel *p_cell_cycle_model =
					new ChondrocytesOnlyCellCycleModel;
				double birth_time =
					-p_cell_cycle_model->GetAverageTransitCellCycleTime() * RandomNumberGenerator::Instance()->ranf();
				cell->SetBirthTime(birth_time);
			}
			else
			{
				cell->SetBirthTime(0.0); //Average transit cell cycle time is 12.0 with current values
										  
			}
			//increase counter
			i++;
		}
	}
	mpCellPopulation->AddCellWriter<CellAncestorWriter>();
}

void NodeBasedMesenchymalCondensation::SetDimensions(
	unsigned numberOfCellsWide, unsigned numberOfCellsDeep)
{
	mNumberOfNodesPerXDimension = numberOfCellsWide;
	mNumberOfNodesPerYDimension = numberOfCellsDeep;
}

/**
 * Generates random node coordinates for a 3D cell sheet a given number of cells wide and deep.
 * Arrangement of the nodes will be on a Cartesian grid.
 */
void NodeBasedMesenchymalCondensation::GenerateNodesOnCartesianGrid()
{

	mNodes.clear();
	unsigned n_nodes_width = mNumberOfNodesPerXDimension;
	unsigned n_nodes_depth = mNumberOfNodesPerYDimension;
	unsigned n_nodes = n_nodes_width * n_nodes_depth;
	mNodes.reserve(n_nodes);

	unsigned id = 0;

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

                double z_offset = mDistanceBetweeenBoundaries *RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = i + noise * cos(random_azimuth_angle) * sin(random_zenith_angle);
				double y_coordinate = j + noise * sin(random_azimuth_angle) * sin(random_zenith_angle);
				double z_coordinate = z_offset + noise * cos(random_zenith_angle);
				mNodes.push_back(
					new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
		}
	}
	
	mNodesGenerated = true;
}

/**
 * Generates random node coordinates for a 3D cell sheet a given number of cells wide, deep and high.
 * Arrangement of the nodes will be on a hcp lattice.
 */
void NodeBasedMesenchymalCondensation::GenerateNodesOnHCPGrid()
{

	mNodes.clear();
	unsigned n_nodes_width = mNumberOfNodesPerXDimension;
	unsigned n_nodes_depth = mNumberOfNodesPerYDimension;
	unsigned n_nodes = n_nodes_width * n_nodes_depth;
	mNodes.reserve(n_nodes);

	unsigned id = 0;
	unsigned k = 0;

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

				double z_offset = mDistanceBetweeenBoundaries *RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				double x_coordinate = (2 * i + ((j + k) % 2)) * 0.5 + noise * cos(random_azimuth_angle) * sin(random_zenith_angle);
				double y_coordinate = (sqrt(3) * (j + (k % 2) / 3.0)) * 0.5 + noise * sin(random_azimuth_angle) * sin(random_zenith_angle);
				double z_coordinate = z_offset + (2 * sqrt(6) * k / 3.0) * 0.5 + noise * cos(random_zenith_angle);
				mNodes.push_back(
					new Node<3>(id, false, x_coordinate, y_coordinate,
								z_coordinate));
				id++;
			}
	}
	mNodesGenerated = true;
}


double NodeBasedMesenchymalCondensation::getMaxCoordinatePerturbation() const
{
	return mMaxCoordinatePerturbation;
}

void NodeBasedMesenchymalCondensation::setMaxCoordinatePerturbation(
	double maxCoordinatePerturbation)
{
	mMaxCoordinatePerturbation = maxCoordinatePerturbation;
}

double NodeBasedMesenchymalCondensation::getDistanceBetweeenBoundaries() const
{
	return mDistanceBetweeenBoundaries;
}

void NodeBasedMesenchymalCondensation::setDistanceBetweeenBoundaries(
	double distance)
{
	mDistanceBetweeenBoundaries = distance;
}

unsigned NodeBasedMesenchymalCondensation::getNumberOfNodesPerXDimension() const
{
	return mNumberOfNodesPerXDimension;
}

unsigned NodeBasedMesenchymalCondensation::getNumberOfNodesPerYDimension() const
{
	return mNumberOfNodesPerYDimension;
}

unsigned NodeBasedMesenchymalCondensation::getSeed() const
{
	return mSeed;
}

void NodeBasedMesenchymalCondensation::setSynchronizeCellCycles(
	bool synchronizeCellCycles)
{
	mSynchronizeCellCycles = synchronizeCellCycles;
}

void NodeBasedMesenchymalCondensation::setDivisionDirections(
	bool divisionDirections)
{
	mDivisionDirections = divisionDirections;
}

void NodeBasedMesenchymalCondensation::UseRandomSeed()
{

	// Reseed the number generator so that different runs will actually produce different results
	mSeed = time(NULL);
	RandomNumberGenerator::Instance()->Reseed(mSeed);
}
