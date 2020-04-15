/*
 * AbstractPerichondrialLayer.cpp
 *
 *  Created on: May 19, 2017
 *      Author: Sonja Mathias
 */

#include "AbstractPerichondrialLayer.hpp"
#include "Exception.hpp"

AbstractPerichondrialLayer::AbstractPerichondrialLayer() {
	// Subclasses should always call the other constructor.
	NEVER_REACHED;

	}

AbstractPerichondrialLayer::AbstractPerichondrialLayer(unsigned colour) :
		AbstractCellProperty(), mColour(colour) {

}

AbstractPerichondrialLayer::~AbstractPerichondrialLayer() {
}

unsigned AbstractPerichondrialLayer::GetColour() const
{
    return mColour;
}
