/*
 * FixedCellDivisionDirection.cpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#include "FixedCellDivisionDirection.hpp"

template<unsigned SPACE_DIM> FixedCellDivisionDirection<SPACE_DIM>::FixedCellDivisionDirection() :
		AbstractCellDivisionDirection<SPACE_DIM>(20,
				zero_vector<double>(SPACE_DIM)) {

    if(SPACE_DIM == 3){
        this->mDirection(0) = 0.5;
        this->mDirection(1) = sqrt(3.0)/2.0;
    }
    else{
        this->mDirection(0) = 1.0;
    }	

}

template class FixedCellDivisionDirection<1> ;
template class FixedCellDivisionDirection<2> ;
template class FixedCellDivisionDirection<3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedCellDivisionDirection)
