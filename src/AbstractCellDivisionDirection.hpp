/*
 * AbstractCellDivisionDirection.hpp
 *
 *  Created on: May 15, 2017
 *      Author: Sonja Mathias
 */

#ifndef ABSTRACTCELLDIVISIONDIRECTION_HPP_
#define ABSTRACTCELLDIVISIONDIRECTION_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "AbstractCellProperty.hpp"

#include "UblasVectorInclude.hpp" //needed for c_vector reference


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM> class AbstractCellDivisionDirection: public AbstractCellProperty {

	// Allow tests to access private members to test private functions
	friend class TestAbstractCellDivisionDirection;

private:
	unsigned mColour;


	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
	    archive & boost::serialization::base_object<AbstractCellProperty>(*this);
	    archive & mColour;
	    archive & mDirection;
	}

protected:
	c_vector<double, SPACE_DIM> mDirection;

public:
	AbstractCellDivisionDirection(unsigned colour, c_vector<double, SPACE_DIM> direction);
	virtual ~AbstractCellDivisionDirection();

	unsigned GetColour() const;

	c_vector<double, SPACE_DIM> GetCellDivisionDirection() const;

};

#endif /* ABSTRACTCELLDIVISIONDIRECTION_HPP_ */
