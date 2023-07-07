#ifndef PERTURBEDUPWARDSCELLDIVISIONDIRECTION_HPP_
#define PERTURBEDUPWARDSCELLDIVISIONDIRECTION_HPP_

#include "AbstractCellDivisionDirection.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned SPACE_DIM>
class PerturbedUpwardsCellDivisionDirection: public AbstractCellDivisionDirection<
		SPACE_DIM> {

private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive
				& boost::serialization::base_object<
						AbstractCellDivisionDirection<SPACE_DIM> >(*this);
	}
	
	double mMaximumZenithAngle;
	
public:
	PerturbedUpwardsCellDivisionDirection();
    
    void setMaximumZenithAngle(double);
    double getMaximumZenithAngle() const;
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PerturbedUpwardsCellDivisionDirection)

#endif /* PERTURBEDUPWARDSCELLDIVISIONDIRECTION_HPP_ */
