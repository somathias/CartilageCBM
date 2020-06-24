/*
 * DirectionalRestLengthPWQForce.hpp
 *
 *  Created on: June 24, 2020
 *      Author: Sonja Mathias
 */

#ifndef DIRECTIONALRESTLENGTHPWQFORCE_HPP_
#define DIRECTIONALRESTLENGTHPWQFORCE_HPP_

#include "PWQGeneralisedLinearSpringForce.hpp"


class DirectionalRestLengthPWQForce: public PWQGeneralisedLinearSpringForce<3> {

private:

	/** Needed for serialization. */
	friend class boost::serialization::access;
	/**
	 * Archive the object and its member variables.
	 *
	 * @param archive the archive
	 * @param version the current version of this class
	 */
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive
				& boost::serialization::base_object<
						PWQGeneralisedLinearSpringForce<3> >(
						*this);
	}


public:
	DirectionalRestLengthPWQForce();

	/**
	 * Overridden CalculateForceBetweenNodes method which uses different spring
	 * stiffness values for the repulsion and adhesion parts of the force law.
	 * It also uses the member variable mAlpha instead of a hard-coded decay
	 * in the force function.
	 */
	c_vector<double, 3> CalculateForceBetweenNodes(
			unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
			AbstractCellPopulation<3>& rCellPopulation);

	/**
	 * Overridden OutputForceParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalRestLengthPWQForce)
CHASTE_CLASS_EXPORT(DirectionalRestLengthPWQForce)

#endif /* DIRECTIONALRESTLENGTHPWQFORCE_HPP_ */
