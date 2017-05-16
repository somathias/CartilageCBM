/*
 * AbstractCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "AbstractCellDivisionDirection.hpp"

template<unsigned SPACE_DIM> AbstractCellDivisionDirection<SPACE_DIM>::AbstractCellDivisionDirection(
		unsigned colour, c_vector<double, SPACE_DIM> direction) :
		AbstractCellProperty(), mColour(colour), mDirection(direction) {

}

template<unsigned SPACE_DIM> AbstractCellDivisionDirection<SPACE_DIM>::~AbstractCellDivisionDirection() {

}

template<unsigned SPACE_DIM>
unsigned AbstractCellDivisionDirection<SPACE_DIM>::GetColour() const {
	return mColour;
}

template<unsigned SPACE_DIM> c_vector<double, SPACE_DIM> AbstractCellDivisionDirection<
		SPACE_DIM>::GetCellDivisionDirection() const {
	return mDirection;
}

template class AbstractCellDivisionDirection<1> ;
template class AbstractCellDivisionDirection<2> ;
template class AbstractCellDivisionDirection<3> ;
