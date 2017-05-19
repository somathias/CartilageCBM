/*
 * AbstractCellTissueType.hpp
 *
 *  Created on: May 19, 2017
 *      Author: Sonja Mathias
 */

#ifndef ABSTRACTCELLTISSUETYPE_HPP_
#define ABSTRACTCELLTISSUETYPE_HPP_

class AbstractCellTissueType: public AbstractCellProperty {

private:
	unsigned mColour;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive
				& boost::serialization::base_object<AbstractCellProperty>(
						*this);
		archive & mColour;
	}

	/**
	 * Default constructor needs to be defined for archiving, but never actually used,
	 * since subclasses call the normal constructor.
	 */
	AbstractCellTissueType();

public:
	/**
	 * Constructor.
	 *
	 * @param colour  what colour cells with this proliferative type should be in the visualizer
	 */
	AbstractCellTissueType(unsigned colour);
	virtual ~AbstractCellTissueType();

	/**
	 * @return #mColour.
	 */
	unsigned GetColour() const;
};

#endif /* ABSTRACTCELLTISSUETYPE_HPP_ */
