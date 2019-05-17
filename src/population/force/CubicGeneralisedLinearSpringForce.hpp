/*
 * CubicGeneralisedLinearSpringForce.hpp
 *
 *  Created on: Apr 30, 2019
 *      Author: Sonja Mathias
 */

#ifndef CUBICGENERALISEDLINEARSPRINGFORCE_HPP_
#define CUBICGENERALISEDLINEARSPRINGFORCE_HPP_

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A force class extending the GeneralisedLinearSpringForce class
 * with the functionality of setting the spring stiffness values
 * for the repulsion and adhesion parts individually.
 *
 * mMeinekeSpringStiffness -> adhesion part (ie no overlap between cells)
 * mRepulsionSpringStiffness -> repulsion part (cells overlap)
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class CubicGeneralisedLinearSpringForce: public GeneralisedLinearSpringForce<
		ELEMENT_DIM, SPACE_DIM> {

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
						GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM> >(
						*this);
	}


public:
	CubicGeneralisedLinearSpringForce();


	/**
	 * Overridden CalculateForceBetweenNodes method which uses different spring
	 * stiffness values for the repulsion and adhesion parts of the force law.
	 * It also uses the member variable mAlpha instead of a hard-coded decay
	 * in the force function.
	 */
	c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(
			unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
			AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

	/**
	 * Overridden OutputForceParameters() method.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(
		CubicGeneralisedLinearSpringForce)

#endif /* CUBICGENERALISEDLINEARSPRINGFORCE_HPP_ */
