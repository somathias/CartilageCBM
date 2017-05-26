/*
 * ChondrocyteCellTissueType.cpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#include "ChondrocyteCellTissueType.hpp"

ChondrocyteCellTissueType::ChondrocyteCellTissueType() :
		AbstractCellTissueType(24) {

}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ChondrocyteCellTissueType)
