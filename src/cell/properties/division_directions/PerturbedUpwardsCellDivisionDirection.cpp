#include "PerturbedUpwardsCellDivisionDirection.hpp"

template<unsigned SPACE_DIM> PerturbedUpwardsCellDivisionDirection<SPACE_DIM>::PerturbedUpwardsCellDivisionDirection() :
		AbstractCellDivisionDirection<SPACE_DIM>(21,
				zero_vector<double>(SPACE_DIM)), mMaximumZenithAngle(0.25*M_PI){
//	c_vector<double, SPACE_DIM> main_direction = zero_vector<double>(SPACE_DIM);
//	main_direction(SPACE_DIM - 1) = 1.0;          //upwards
//	AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM>(17, main_direction);

	this->mDirection(SPACE_DIM - 1) = 1.0;

}


template<unsigned SPACE_DIM> void PerturbedUpwardsCellDivisionDirection<SPACE_DIM>::setMaximumZenithAngle(double theta_max) 
{
    //theta_max should be between zero and pi. 
    assert(theta_max > 0);
    assert(theta_max < M_PI);
    mMaximumZenithAngle = theta_max;
}

template<unsigned SPACE_DIM> double PerturbedUpwardsCellDivisionDirection<SPACE_DIM>::getMaximumZenithAngle() const 
{
    return mMaximumZenithAngle;
}


template class PerturbedUpwardsCellDivisionDirection<1> ;
template class PerturbedUpwardsCellDivisionDirection<2> ;
template class PerturbedUpwardsCellDivisionDirection<3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PerturbedUpwardsCellDivisionDirection)
