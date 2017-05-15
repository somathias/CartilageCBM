/*
 * UpwardsCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "UpwardsCellDivisionDirection.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> UpwardsCellDivisionDirection<
		ELEMENT_DIM, SPACE_DIM>::UpwardsCellDivisionDirection() :
		AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>(17,
				zero_vector<double>(SPACE_DIM)) {
//	c_vector<double, SPACE_DIM> main_direction = zero_vector<double>(SPACE_DIM);
//	main_direction(SPACE_DIM - 1) = 1.0;          //upwards
//	AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>(17, main_direction);

	this->mDirection(SPACE_DIM - 1) = 1.0;

}

template class UpwardsCellDivisionDirection<1,1>;
template class UpwardsCellDivisionDirection<1,2>;
template class UpwardsCellDivisionDirection<2,2>;
template class UpwardsCellDivisionDirection<1,3>;
template class UpwardsCellDivisionDirection<2,3>;
template class UpwardsCellDivisionDirection<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(UpwardsCellDivisionDirection)
