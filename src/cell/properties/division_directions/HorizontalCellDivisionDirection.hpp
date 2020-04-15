/*
 * HorizontalCellDivisionDirection.hpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#ifndef HORIZONTALCELLDIVISIONDIRECTION_HPP_
#define HORIZONTALCELLDIVISIONDIRECTION_HPP_

#include "AbstractCellDivisionDirection.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned SPACE_DIM>
class HorizontalCellDivisionDirection: public AbstractCellDivisionDirection<
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
	HorizontalCellDivisionDirection();
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HorizontalCellDivisionDirection)

#endif /*HORIZONTALCELLDIVISIONDIRECTION_HPP_ */

