/*
 * DirectionalAdhesionForce.hpp
 *
 *  Created on: Feb 26, 2018
 *      Author: Sonja Mathias
 */

#ifndef DIRECTIONALADHESIONFORCE_HPP_
#define DIRECTIONALADHESIONFORCE_HPP_

#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class DirectionalAdhesionForce: public CellTissueTypeBasedGeneralisedLinearSpringForce<
		ELEMENT_DIM, SPACE_DIM> {

private:

	/**
	 * Baseline adhesion multiplier.
	 *
	 * Defaults to 0.1 in the constructor.
	 */
	double mBaselineAdhesionMultiplier;

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
						CellTissueTypeBasedGeneralisedLinearSpringForce<
								ELEMENT_DIM, SPACE_DIM> >(*this);
		archive & mBaselineAdhesionMultiplier;
	}

public:
	DirectionalAdhesionForce();

	/**
	 * Overridden VariableSpringConstantMultiplicationFactor() method.
	 *
	 * This method takes account of the distinct spring constants present
	 * for homotypic (perichondrial-perichondrial and chondrocyte-chondrocyte) and
	 * heterotypic (perichondrial-chondrocyte) interactions between neighbouring
	 * cells. Additionally, it also scales the spring constants by the directional
	 * alignment of the interactions with the main adhesion direction for the
	 * different cell tissue types. This means that homotypic adhesion between
	 * perichondrial cells is strongest in the (x,y) plane and homotypic adhesion
	 * between chondrocytes is strongest perpendicular to the (x,y) plane.
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
			AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
			bool isCloserThanRestLength);
	/**
	 * Overridden OutputForceParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalAdhesionForce)

#endif /* DIRECTIONALADHESIONFORCE_HPP_ */
