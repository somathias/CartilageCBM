/*
 * FixedCellDivisionDirection.hpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#ifndef FIXEDCELLDIVISIONDIRECTION_HPP_
#define FIXEDCELLDIVISIONDIRECTION_HPP_

#include "AbstractCellDivisionDirection.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned SPACE_DIM>
class FixedCellDivisionDirection: public AbstractCellDivisionDirection<
		SPACE_DIM> {

private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive
				& boost::serialization::base_object<
						AbstractCellDivisionDirection<SPACE_DIM> >(*this);
	}
public:
	FixedCellDivisionDirection();
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedCellDivisionDirection)

#endif /*FIXEDCELLDIVISIONDIRECTION_HPP_ */

