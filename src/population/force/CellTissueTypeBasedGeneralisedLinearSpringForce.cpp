/*
 * CellTissueTypeBasedGeneralisedLinearSpringForce.cpp
 *
 *  Created on: May 30, 2017
 *      Author: kubuntu1404
 */

#include "CellTissueTypeBasedGeneralisedLinearSpringForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::CellTissueTypeBasedGeneralisedLinearSpringForce() :
		GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(), mHomotypicPerichondrialSpringConstantMultiplier(
				1.0), mHomotypicChondrocyteSpringConstantMultiplier(1.0), mHeterotypicSpringConstantMultiplier(
				1.0) {
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
		unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
		bool isCloserThanRestLength) {

	if (isCloserThanRestLength) {
		return 1.0;
	} else {
		// Determine which type the cells corresponding to these nodes are
		CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(
				nodeAGlobalIndex);
		bool cell_A_is_perichondrial = p_cell_A->template HasCellProperty<
				PerichondrialCellTissueType>();
		bool cell_A_is_chondrocyte = p_cell_A->template HasCellProperty<
				ChondrocyteCellTissueType>();

		CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(
				nodeBGlobalIndex);
		bool cell_B_is_perichondrial = p_cell_B->template HasCellProperty<
				PerichondrialCellTissueType>();
		bool cell_B_is_chondrocyte = p_cell_B->template HasCellProperty<
				ChondrocyteCellTissueType>();


		if ((cell_A_is_perichondrial || cell_A_is_chondrocyte)
				&& (cell_B_is_perichondrial || cell_B_is_chondrocyte)) {
			//both of them have a cell tissue type

			if (cell_A_is_perichondrial && cell_B_is_perichondrial) {
				// For homotypic interactions between perichondrial cells, scale the spring constant by mHomotypicPerichondrialSpringConstantMultiplier
				return mHomotypicPerichondrialSpringConstantMultiplier;
			} else if (cell_A_is_chondrocyte && cell_B_is_chondrocyte) {
				// For homotypic interactions between chondrocyte cells, scale the spring constant by mHomotypicChondrocyteSpringConstantMultiplier
				return mHomotypicChondrocyteSpringConstantMultiplier;
			} else {
				// For heterotypic interactions, scale the spring constant by mHeterotypicSpringConstantMultiplier
				return mHeterotypicSpringConstantMultiplier;
			}
		} else {
			//leave the spring constant unchanged if one of them does not have a cell tissue type
			return 1.0;
		}
	}
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHomotypicPerichondrialSpringConstantMultiplier() {
	return mHomotypicPerichondrialSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHomotypicPerichondrialSpringConstantMultiplier(
		double perichondrialSpringConstantMultiplier) {
	assert(perichondrialSpringConstantMultiplier > 0.0);
	mHomotypicPerichondrialSpringConstantMultiplier =
			perichondrialSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHomotypicChondrocyteSpringConstantMultiplier() {
	return mHomotypicChondrocyteSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHomotypicChondrocyteSpringConstantMultiplier(
		double chondrocyteSpringConstantMultiplier) {
	assert(chondrocyteSpringConstantMultiplier > 0.0);
	mHomotypicChondrocyteSpringConstantMultiplier =
			chondrocyteSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetHeterotypicSpringConstantMultiplier() {
	return mHeterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetHeterotypicSpringConstantMultiplier(
		double heterotypicSpringConstantMultiplier) {
	assert(heterotypicSpringConstantMultiplier > 0.0);
	mHeterotypicSpringConstantMultiplier = heterotypicSpringConstantMultiplier;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(
		out_stream& rParamsFile) {
	*rParamsFile << "\t\t\t<HomotypicPerichondrialSpringConstantMultiplier>"
			<< mHomotypicPerichondrialSpringConstantMultiplier
			<< "</HomotypicPerichondrialSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<HomotypicChondrocyteSpringConstantMultiplier>"
			<< mHomotypicChondrocyteSpringConstantMultiplier
			<< "</HomotypicChondrocyteSpringConstantMultiplier>\n";
	*rParamsFile << "\t\t\t<HeterotypicSpringConstantMultiplier>"
			<< mHeterotypicSpringConstantMultiplier
			<< "</HeterotypicSpringConstantMultiplier>\n";

	// Call direct parent class
	GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(
			rParamsFile);
}

// Explicit instantiation
template class CellTissueTypeBasedGeneralisedLinearSpringForce<1, 1> ;
template class CellTissueTypeBasedGeneralisedLinearSpringForce<1, 2> ;
template class CellTissueTypeBasedGeneralisedLinearSpringForce<2, 2> ;
template class CellTissueTypeBasedGeneralisedLinearSpringForce<1, 3> ;
template class CellTissueTypeBasedGeneralisedLinearSpringForce<2, 3> ;
template class CellTissueTypeBasedGeneralisedLinearSpringForce<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellTissueTypeBasedGeneralisedLinearSpringForce)

