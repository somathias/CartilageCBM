/*
 * AbstractCellTissueType.cpp
 *
 *  Created on: May 19, 2017
 *      Author: Sonja Mathias
 */

#include "AbstractCellTissueType.hpp"
#include "Exception.hpp"

AbstractCellTissueType::AbstractCellTissueType() {
	// Subclasses should always call the other constructor.
	NEVER_REACHED;

	}

AbstractCellTissueType::AbstractCellTissueType(unsigned colour) :
		AbstractCellProperty(), mColour(colour) {

}

AbstractCellTissueType::~AbstractCellTissueType() {
}

unsigned AbstractCellTissueType::GetColour() const
{
    return mColour;
}
