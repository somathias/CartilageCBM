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
#include <cxxtest/TestSuite.h>
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
#include "CellAncestorWriter.hpp"

#include "CellTissueTypeBasedCellCycleModel.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"

#include "FakePetscSetup.hpp"

class NodeBasedCartilageSheet {

private:
	//c_vector<unsigned, 3> mNumberOfNodesPerDimension;
	unsigned mNumberOfNodesPerXDimension;
	unsigned mNumberOfNodesPerYDimension;
	unsigned mNumberOfNodesPerZDimension;
	double mMaxCoordinatePerturbation;

	unsigned mSeed;

	bool mSynchronizeCellCycles;

	std::vector<Node<3>*> mNodes;
	bool mNodesGenerated;
	NodesOnlyMesh<3> mMesh;
	std::vector<CellPtr> mCells;

	bool mCellPopulationSetup;
	boost::shared_ptr<NodeBasedCellPopulation<3> > mpCellPopulation;

public:

	NodeBasedCartilageSheet();
	virtual ~NodeBasedCartilageSheet();

	void Setup() throw (Exception);
	bool isCellPopulationSetup() const;
	boost::shared_ptr<NodeBasedCellPopulation<3> > GetCellPopulation()
			throw (Exception);

	void InitialiseTissueLayersAndCellDivisionDirections() throw (Exception);
	void InitialiseBulkStemCellConfiguration(unsigned numberOfCellsWide,
			unsigned numberOfCellsDeep) throw (Exception);

	void SetCartilageSheetDimensions(unsigned, unsigned, unsigned);

	void GenerateNodesOnCartesianGrid();
	void GenerateNodesOnHCPGrid();

	void UseRandomSeed();
	double getMaxCoordinatePerturbation() const;
	void setMaxCoordinatePerturbation(double maxCoordinatePerturbation);
	unsigned getNumberOfNodesPerXDimension() const;
	unsigned getNumberOfNodesPerYDimension() const;
	unsigned getNumberOfNodesPerZDimension() const;
	unsigned getSeed() const;
	void setSynchronizeCellCycles(bool synchronizeCellCycles);
};

#endif /* NODEBASEDCARTILAGESHEET_HPP_ */
