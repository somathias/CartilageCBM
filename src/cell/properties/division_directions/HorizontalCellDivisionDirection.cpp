/*
 * HorizontalCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "HorizontalCellDivisionDirection.hpp"

template<unsigned SPACE_DIM> HorizontalCellDivisionDirection<SPACE_DIM>::HorizontalCellDivisionDirection() :
		AbstractCellDivisionDirection<SPACE_DIM>(19,
				zero_vector<double>(SPACE_DIM)) {
}

template class HorizontalCellDivisionDirection<1> ;
template class HorizontalCellDivisionDirection<2> ;
template class HorizontalCellDivisionDirection<3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HorizontalCellDivisionDirection)