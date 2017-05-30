/*
 * PerichondrialCellTissueType.cpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#include "PerichondrialCellTissueType.hpp"

PerichondrialCellTissueType::PerichondrialCellTissueType() :
		AbstractCellTissueType(23) {

}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PerichondrialCellTissueType)
