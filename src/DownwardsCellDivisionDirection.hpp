/*
 * DownwardsCellDivisionDirection.hpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#ifndef DOWNWARDSCELLDIVISIONDIRECTION_HPP_
#define DOWNWARDSCELLDIVISIONDIRECTION_HPP_

#include "AbstractCellDivisionDirection.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class DownwardsCellDivisionDirection: public AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM> {

private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive	& boost::serialization::base_object<AbstractCellDivisionDirection<ELEMENT_DIM, SPACE_DIM> > (*this);
	}
public:
	DownwardsCellDivisionDirection();
};




#endif /* DOWNWARDSCELLDIVISIONDIRECTION_HPP_ */
