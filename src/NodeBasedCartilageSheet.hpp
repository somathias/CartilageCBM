/*
 * NodeBasedCartilageSheet.hpp
 *
 *  Created on: Jul 25, 2017
 *      Author: kubuntu1404
 */

#ifndef NODEBASEDCARTILAGESHEET_HPP_
#define NODEBASEDCARTILAGESHEET_HPP_

#include <ctime>
#include <string>
#include <sstream>
#include <iostream>
#include "CheckpointArchiveTypes.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellAncestor.hpp"
#include "CellAncestorWriter.hpp"
#include "WildTypeCellMutationState.hpp"


#include "CellTissueTypeBasedCellCycleModel.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"
#include "OrientationBasedDivisionRule.hpp"


class NodeBasedCartilageSheet {

	// Allow tests to access private members to test private functions
	friend class TestNodeBasedCartilageSheet;

private:
	//c_vector<unsigned, 3> mNumberOfNodesPerDimension;
	unsigned mNumberOfNodesPerXDimension;
	unsigned mNumberOfNodesPerYDimension;
	unsigned mNumberOfNodesPerZDimension;

	unsigned mNumberOfPerichondrialLayersAbove;
	unsigned mNumberOfPerichondrialLayersBelow;

	double mMaxCoordinatePerturbation;

	unsigned mSeed;

	bool mSynchronizeCellCycles;
	bool mDivisionDirections;

	std::vector<Node<3>*> mNodes;
	bool mNodesGenerated;
	NodesOnlyMesh<3> mMesh;
	std::vector<CellPtr> mCells;

	bool mCellPopulationSetup;
	boost::shared_ptr<NodeBasedCellPopulation<3> > mpCellPopulation;

public:

	NodeBasedCartilageSheet();
	virtual ~NodeBasedCartilageSheet();

	void Setup();
	bool isCellPopulationSetup() const;
	boost::shared_ptr<NodeBasedCellPopulation<3> > GetCellPopulation();

	void InitialiseTissueLayersAndCellDivisionDirections();
	void InitialiseBulkStemCellConfiguration(unsigned, unsigned);
	void InitialiseRandomStemCellConfiguration(unsigned);

	void SetCartilageSheetDimensions(unsigned, unsigned, unsigned);

	void GenerateNodesOnCartesianGrid();
	void GenerateNodesOnHCPGrid();

	void UseRandomSeed();

	double getMaxCoordinatePerturbation() const;
	void setMaxCoordinatePerturbation(double maxCoordinatePerturbation);

	unsigned getNumberOfNodesPerXDimension() const;
	unsigned getNumberOfNodesPerYDimension() const;
	unsigned getNumberOfNodesPerZDimension() const;

	void setNumberOfPerichondrialLayersAbove(unsigned);
	void setNumberOfPerichondrialLayersBelow(unsigned);

	unsigned getNumberOfPerichondrialLayersAbove() const;
	unsigned getNumberOfPerichondrialLayersBelow() const;

	unsigned getSeed() const;

	void setSynchronizeCellCycles(bool synchronizeCellCycles);
	void setDivisionDirections(bool divisionDirections);
};

#endif /* NODEBASEDCARTILAGESHEET_HPP_ */
