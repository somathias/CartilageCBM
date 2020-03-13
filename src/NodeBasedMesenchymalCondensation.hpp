/*
 * NodeBasedMesenchymalCondensation.hpp
 *
 *  Created on: Mar 12, 2020
 *      Author: Sonja Mathias
 */

#ifndef NODEBASEDMESENCHYMALCONDENSATION_HPP_
#define NODEBASEDMESENCHYMALCONDENSATION_HPP_

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


#include "ChondrocytesOnlyCellCycleModel.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"
#include "OrientationBasedDivisionRule.hpp"


class NodeBasedMesenchymalCondensation {

	// Allow tests to access private members to test private functions
	friend class TestNodeBasedMesenchymalCondensation;

private:
	unsigned mNumberOfNodesPerXDimension;
	unsigned mNumberOfNodesPerYDimension;

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

	NodeBasedMesenchymalCondensation();
	virtual ~NodeBasedMesenchymalCondensation();

	void Setup();
	bool isCellPopulationSetup() const;
	boost::shared_ptr<NodeBasedCellPopulation<3> > GetCellPopulation();

	void InitialiseRandomConfiguration(unsigned);

	void SetDimensions(unsigned, unsigned);

	void GenerateNodesOnCartesianGrid();
	void GenerateNodesOnHCPGrid();

	void UseRandomSeed();

	double getMaxCoordinatePerturbation() const;
	void setMaxCoordinatePerturbation(double maxCoordinatePerturbation);

	unsigned getNumberOfNodesPerXDimension() const;
	unsigned getNumberOfNodesPerYDimension() const;

	unsigned getSeed() const;

	void setSynchronizeCellCycles(bool synchronizeCellCycles);
	void setDivisionDirections(bool divisionDirections);
};

#endif /* NODEBASEDMESENCHYMALCONDENSATION_HPP_ */
