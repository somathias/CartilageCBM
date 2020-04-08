/*
 * AbstractPerichondrialLayer.hpp
 *
 *  Created on: May 19, 2017
 *      Author: Sonja Mathias
 */

#ifndef ABSTRACTPERICHONDRIALLAYER_HPP_
#define ABSTRACTPERICHONDRIALLAYER_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class AbstractPerichondrialLayer : public AbstractCellProperty {

private:
	unsigned mColour;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive	& boost::serialization::base_object<AbstractCellProperty>(
						*this);
		archive & mColour;
	}

	/**
	 * Default constructor needs to be defined for archiving, but never actually used,
	 * since subclasses call the normal constructor.
	 */
	AbstractPerichondrialLayer();

public:


	/**
	 * Constructor.
	 *
	 * @param colour  what colour cells with this proliferative type should be in the visualizer
	 */
	AbstractPerichondrialLayer(unsigned colour);
	virtual ~AbstractPerichondrialLayer();

	/**
	 * @return #mColour.
	 */
	unsigned GetColour() const;
};

#endif /* ABSTRACTPERICHONDRIALLAYER_HPP_ */
