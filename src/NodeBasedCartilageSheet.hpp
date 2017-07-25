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
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellAncestorWriter.hpp"
#include "StochasticDurationCellCycleModel.hpp"

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
	NodesOnlyMesh<3> mesh;
	std::vector< CellPtr > cells;

public:

	NodeBasedCartilageSheet();
	virtual ~NodeBasedCartilageSheet();

	boost::shared_ptr<NodeBasedCellPopulation<3> > Setup();

	void GenerateRandomHCPNodes(
			std::vector<Node<3>*> & rNodes, unsigned n_nodes_width,
			unsigned n_nodes_depth, unsigned n_nodes_height, double max_noise);

};

#endif /* NODEBASEDCARTILAGESHEET_HPP_ */
