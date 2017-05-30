/*
 * CellTissueBasedGeneralisedLinearSpringForce.hpp
 *
 *  Created on: May 30, 2017
 *      Author: Sonja Mathias
 */

#ifndef CELLTISSUEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_
#define CELLTISSUEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A force class copying the DifferentialAdhesionGeneralisedSpringForce class,
 * but based on cell tissue types instead of cell labels.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class CellTissueTypeBasedGeneralisedLinearSpringForce: public GeneralisedLinearSpringForce<
		ELEMENT_DIM, SPACE_DIM> {

private:

	/**
	 * A scalar determining the relative spring constant for homotypic
	 * interactions between neighbouring perichondrial cells, used in the
	 * overridden method VariableSpringConstantMultiplicationFactor().
	 *
	 * Defaults to 1.0 in the constructor.
	 */
	double mHomotypicPerichondrialSpringConstantMultiplier;

	/**
	 * A scalar determining the relative spring constant for homotypic
	 * interactions between neighbouring chondrocyte cells, used in the
	 * overridden method VariableSpringConstantMultiplicationFactor().
	 *
	 * Defaults to 1.0 in the constructor.
	 */
	double mHomotypicChondrocyteSpringConstantMultiplier;

	/**
	 * A scalar determining the relative spring constant for heterotypic
	 * (perichondrial-chondrocyte) interactions between neighbouring cells, used
	 * in the overridden method VariableSpringConstantMultiplicationFactor().
	 *
	 * Defaults to 1.0 in the constructor.
	 */
	double mHeterotypicSpringConstantMultiplier;

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
						GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM> >(
						*this);
		archive & mHomotypicPerichondrialSpringConstantMultiplier;
		archive & mHomotypicChondrocyteSpringConstantMultiplier;
		archive & mHeterotypicSpringConstantMultiplier;
	}
public:
	/**
	 * Constructor.
	 */
	CellTissueTypeBasedGeneralisedLinearSpringForce();

	/**
	 * Overridden VariableSpringConstantMultiplicationFactor() method.
	 *
	 * This method takes account of the distinct spring constants present
	 * for homotypic (perichondrial-perichondrial and chondrocyte-chondrocyte) and
	 * heterotypic (perichondrial-chondrocyte) interactions between neighbouring
	 * cells.
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
	 * @return #mHomotypicPerichondrialSpringConstantMultiplier.
	 */
	double GetHomotypicPerichondrialSpringConstantMultiplier();

	/**
	 * Set mHomotypicPerichondrialSpringConstantMultiplier.
	 *
	 * @param perichondrialSpringConstantMultiplier the new value of mHomotypicPerichondrialSpringConstantMultiplier
	 */
	void SetHomotypicPerichondrialSpringConstantMultiplier(
			double perichondrialSpringConstantMultiplier);

	/**
	 * @return #mHomotypicChondrocyteSpringConstantMultiplier.
	 */
	double GetHomotypicChondrocyteSpringConstantMultiplier();

	/**
	 * Set mHomotypicChondrocyteSpringConstantMultiplier.
	 *
	 * @param chondrocyteSpringConstantMultiplier the new value of mHomotypicChondrocyteSpringConstantMultiplier
	 */
	void SetHomotypicChondrocyteSpringConstantMultiplier(
			double chondrocyteSpringConstantMultiplier);

	/**
	 * @return #mHeterotypicSpringConstantMultiplier.
	 */
	double GetHeterotypicSpringConstantMultiplier();

	/**
	 * Set mHeterotypicSpringConstantMultiplier.
	 *
	 * @param heterotypicSpringConstantMultiplier the new value of mHeterotypicSpringConstantMultiplier
	 */
	void SetHeterotypicSpringConstantMultiplier(
			double heterotypicSpringConstantMultiplier);

	/**
	 * Overridden OutputForceParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellTissueTypeBasedGeneralisedLinearSpringForce)

#endif /* CELLTISSUEBASEDGENERALISEDLINEARSPRINGFORCE_HPP_ */
