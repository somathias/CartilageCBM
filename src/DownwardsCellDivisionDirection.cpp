/*
 * DownwardsCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "DownwardsCellDivisionDirection.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> DownwardsCellDivisionDirection<
		ELEMENT_DIM, SPACE_DIM>::DownwardsCellDivisionDirection() :
		AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>(18,
				zero_vector<double>(SPACE_DIM)) {
//	c_vector<double, SPACE_DIM> main_direction = zero_vector<double>(SPACE_DIM);
//	main_direction(SPACE_DIM - 1) = 1.0;          //upwards
//	AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>(17, main_direction);

	this->mDirection(SPACE_DIM-1) = -1.0;

}

template class DownwardsCellDivisionDirection<1,1>;
template class DownwardsCellDivisionDirection<1,2>;
template class DownwardsCellDivisionDirection<2,2>;
template class DownwardsCellDivisionDirection<1,3>;
template class DownwardsCellDivisionDirection<2,3>;
template class DownwardsCellDivisionDirection<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DownwardsCellDivisionDirection)
