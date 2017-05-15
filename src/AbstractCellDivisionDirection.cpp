/*
 * AbstractCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "AbstractCellDivisionDirection.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> AbstractCellDivisionDirection<
		ELEMENT_DIM, SPACE_DIM>::AbstractCellDivisionDirection(unsigned colour,
		c_vector<double, SPACE_DIM> direction) :
		AbstractCellProperty(), mColour(colour), mDirection(direction) {

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> AbstractCellDivisionDirection<
		ELEMENT_DIM, SPACE_DIM>::~AbstractCellDivisionDirection() {

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>::GetColour() const {
	return mColour;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> c_vector<double, SPACE_DIM> AbstractCellDivisionDirection<
		ELEMENT_DIM, SPACE_DIM>::GetCellDivisionDirection() const {
	return mDirection;
}

template class AbstractCellDivisionDirection<1,1>;
template class AbstractCellDivisionDirection<1,2>;
template class AbstractCellDivisionDirection<2,2>;
template class AbstractCellDivisionDirection<1,3>;
template class AbstractCellDivisionDirection<2,3>;
template class AbstractCellDivisionDirection<3,3>;
