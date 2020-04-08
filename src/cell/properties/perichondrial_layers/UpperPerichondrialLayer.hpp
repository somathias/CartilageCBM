/*
 * UpperPerichondrialLayer.hpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#ifndef UPPERPERICHONDRIALLAYER_HPP_
#define UPPERPERICHONDRIALLAYER_HPP_

#include "AbstractCellTissueType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class UpperPerichondrialLayer: public AbstractCellTissueType {
private:
	/** Needed for serialization. */
	friend class boost::serialization::access;
	/**
	 * Archive the cell tissue type.
	 *
	 * @param archive the archive
	 * @param version the current version of this class
	 */
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive
				& boost::serialization::base_object<
						AbstractCellTissueType>(*this);
	}

public:
	UpperPerichondrialLayer();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(UpperPerichondrialLayer)


#endif /* UPPERPERICHONDRIALLAYER_HPP_ */
