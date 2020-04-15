/*
 * UpperPerichondrialLayer.cpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#include "UpperPerichondrialLayer.hpp"

UpperPerichondrialLayer::UpperPerichondrialLayer() :
		AbstractCellTissueType(46) {

}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(UpperPerichondrialLayer)
