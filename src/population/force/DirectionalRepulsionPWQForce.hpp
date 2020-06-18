/*
 * DirectionalRepulsionPWQForce.hpp
 *
 *  Created on: June 10, 2020
 *      Author: Sonja Mathias
 */

#ifndef DIRECTIONALREPULSIONPWQFORCE_HPP_
#define DIRECTIONALREPULSIONPWQFORCE_HPP_

#include "PWQGeneralisedLinearSpringForce.hpp"

class DirectionalRepulsionPWQForce: public PWQGeneralisedLinearSpringForce<3> {

private:

	/**
	 * Baseline repulsion multiplier.
	 *
	 * Defaults to 0.1 in the constructor.
	 */
	double mBaselineRepulsionMultiplier;

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
						PWQGeneralisedLinearSpringForce<3> >(*this);
		archive & mBaselineRepulsionMultiplier;
	}

public:
	DirectionalRepulsionPWQForce();

	/**
	 * @return #mBaselineRepulsionMultiplier.
	 */
	double GetBaselineRepulsionMultiplier();

	/**
	 * Set mBaselineRepulsionMultiplier.
	 *
	 * @param baselineRepulsionMultiplier the new value of mBaselineRepulsionMultiplier
	 */
	void SetBaselineRepulsionMultiplier(
			double baselineRepulsionMultiplier);

	/**
	 * Overridden VariableSpringConstantMultiplicationFactor() method.
	 *
	 *
	 * @param nodeAGlobalIndex index of one neighbouring node
	 * @param nodeBGlobalIndex index of the other neighbouring node
	 * @param rCellPopulation the cell population
	 * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
	 *
	 * @return the multiplication factor.
	 */
	double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
			unsigned nodeBGlobalIndex,
			AbstractCellPopulation<3>& rCellPopulation,
			bool isCloserThanRestLength);
	/**
	 * Overridden OutputForceParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalRepulsionPWQForce)
CHASTE_CLASS_EXPORT(DirectionalRepulsionPWQForce)

#endif /* DIRECTIONALREPULSIONPWQFORCE_HPP_ */
