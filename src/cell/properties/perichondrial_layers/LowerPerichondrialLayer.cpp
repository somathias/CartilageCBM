/*
 * LowerPerichondrialLayer.cpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#include "LowerPerichondrialLayer.hpp"

LowerPerichondrialLayer::LowerPerichondrialLayer() :
		AbstractCellTissueType(45) {

}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(LowerPerichondrialLayer)
