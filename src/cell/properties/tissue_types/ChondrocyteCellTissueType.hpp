/*
 * ChondrocyteCellTissueType.hpp
 *
 *  Created on: May 26, 2017
 *      Author: Sonja Mathias
 */

#ifndef CHONDROCYTECELLTISSUETYPE_HPP_
#define CHONDROCYTECELLTISSUETYPE_HPP_

#include "AbstractCellTissueType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class ChondrocyteCellTissueType: public AbstractCellTissueType {
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
	ChondrocyteCellTissueType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ChondrocyteCellTissueType)


#endif /* CHONDROCYTECELLTISSUETYPE_HPP_ */
