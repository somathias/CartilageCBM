/*
 * DirectionalAdhesionForce.cpp
 *
 *  Created on: Feb 26, 2018
 *      Author: Sonja Mathias
 */

#include "DirectionalAdhesionForce.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DirectionalAdhesionForce<ELEMENT_DIM, SPACE_DIM>::DirectionalAdhesionForce() :
		CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(), mBaselineAdhesionMultiplier(
				0.1) {

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DirectionalAdhesionForce<ELEMENT_DIM, SPACE_DIM>::VariableSpringConstantMultiplicationFactor(
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

		//test that cell A is not both a perichondrial cell and a chondrocyte
		assert(!(cell_A_is_perichondrial && cell_A_is_chondrocyte));

		CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(
				nodeBGlobalIndex);
		bool cell_B_is_perichondrial = p_cell_B->template HasCellProperty<
				PerichondrialCellTissueType>();
		bool cell_B_is_chondrocyte = p_cell_B->template HasCellProperty<
				ChondrocyteCellTissueType>();

		//test that cell B is not both a perichondrial cell and a chondrocyte
		assert(!(cell_B_is_perichondrial && cell_B_is_chondrocyte));

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
void DirectionalAdhesionForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(
		out_stream& rParamsFile) {

	*rParamsFile << "\t\t\t<BaselineAdhesionMultiplier>"
			<< mBaselineAdhesionMultiplier << "</BaselineAdhesionMultiplier>\n";

	// Call direct parent class
	CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(
			rParamsFile);
}

// Explicit instantiation
template class DirectionalAdhesionForce<1, 1> ;
template class DirectionalAdhesionForce<1, 2> ;
template class DirectionalAdhesionForce<2, 2> ;
template class DirectionalAdhesionForce<1, 3> ;
template class DirectionalAdhesionForce<2, 3> ;
template class DirectionalAdhesionForce<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalAdhesionForce)

