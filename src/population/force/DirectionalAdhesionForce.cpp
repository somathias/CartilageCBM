/*
 * DirectionalAdhesionForce.cpp
 *
 *  Created on: Feb 26, 2018
 *      Author: Sonja Mathias
 */

#include "DirectionalAdhesionForce.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "NodeBasedCellPopulation.hpp"

DirectionalAdhesionForce::DirectionalAdhesionForce() :
		CellTissueTypeBasedGeneralisedLinearSpringForce<3>(), mBaselineAdhesionMultiplier(
				0.1) {

}

double DirectionalAdhesionForce::GetBaselineAdhesionMultiplier() {
	return mBaselineAdhesionMultiplier;
}


void DirectionalAdhesionForce::SetBaselineAdhesionMultiplier(
		double baselineAdhesionMultiplier) {
	assert(baselineAdhesionMultiplier >= 0.0);
	mBaselineAdhesionMultiplier =
			baselineAdhesionMultiplier;
}

double DirectionalAdhesionForce::VariableSpringConstantMultiplicationFactor(
		unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
		AbstractCellPopulation<3>& rCellPopulation,
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

			/*
			 * We now calculate the unit difference vector
			 */
			Node<3>* p_node_a = rCellPopulation.GetNode(
					nodeAGlobalIndex);
			Node<3>* p_node_b = rCellPopulation.GetNode(
					nodeBGlobalIndex);

			// Get the node locations
			c_vector<double, 3> node_a_location =
					p_node_a->rGetLocation();
			c_vector<double, 3> node_b_location =
					p_node_b->rGetLocation();

			// Get the node radii for a NodeBasedCellPopulation
			double node_a_radius = 0.0;
			double node_b_radius = 0.0;

			if (bool(
					dynamic_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation))) {
				node_a_radius = p_node_a->GetRadius();
				node_b_radius = p_node_b->GetRadius();
			}

			// Get the unit vector parallel to the line joining the two nodes
			c_vector<double, 3> unit_difference;
			/*
			 * We use the mesh method GetVectorFromAtoB() to compute the direction of the
			 * unit vector along the line joining the two nodes, rather than simply subtract
			 * their positions, because this method can be overloaded (e.g. to enforce a
			 * periodic boundary in Cylindrical2dMesh).
			 */
			unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(
					node_a_location, node_b_location);

			// Calculate the distance between the two nodes
			double distance_between_nodes = norm_2(unit_difference);
			assert(distance_between_nodes > 0);
			assert(!std::isnan(distance_between_nodes));

			unit_difference /= distance_between_nodes;

			/**
			 * Then calculate the dot product between this vector and the unit normal to the x,y plane, n_z = (0,0,1).
			 * Hence this is equal to the z-coordinate of the unit_difference.
			 */
//			c_vector<double, 3> unit_normal_z = zero_vector<double>(3);
//			unit_normal_z[2] = 1.0;

			double abs_dot_product  = abs(unit_difference[2]);

			if (cell_A_is_perichondrial && cell_B_is_perichondrial) {
				/*
				 * For homotypic interactions between perichondrial cells, scale the spring constant
				 * by mHomotypicPerichondrialSpringConstantMultiplier times the directional influence
				 * which should be strongest parallel to the x,y plane (abs_dot_product == 0.0)
				 */
				double directional_multiplier = (1.0-abs_dot_product) * (1-mBaselineAdhesionMultiplier) + mBaselineAdhesionMultiplier;
				return CellTissueTypeBasedGeneralisedLinearSpringForce<3>::mHomotypicPerichondrialSpringConstantMultiplier * directional_multiplier;
			} else if (cell_A_is_chondrocyte && cell_B_is_chondrocyte) {
				/*
				 * For homotypic interactions between chondrocyte cells, scale the spring constant
				 * by mHomotypicChondrocyteSpringConstantMultiplier times the directional influence
				 * which should be strongest perpendicular to the x,y plane (abs_dot_product == 1.0)
				 */
				double directional_multiplier = abs_dot_product * (1-mBaselineAdhesionMultiplier) + mBaselineAdhesionMultiplier;
				return CellTissueTypeBasedGeneralisedLinearSpringForce<3>::mHomotypicChondrocyteSpringConstantMultiplier*directional_multiplier;
			} else {
				// For heterotypic interactions, scale the spring constant by mHeterotypicSpringConstantMultiplier
				return CellTissueTypeBasedGeneralisedLinearSpringForce<3>::mHeterotypicSpringConstantMultiplier;
			}
		} else {
			//leave the spring constant unchanged if one of them does not have a cell tissue type
			return 1.0;
		}
	}
}

void DirectionalAdhesionForce::OutputForceParameters(out_stream& rParamsFile) {

	*rParamsFile << "\t\t\t<BaselineAdhesionMultiplier>"
			<< mBaselineAdhesionMultiplier << "</BaselineAdhesionMultiplier>\n";

	// Call direct parent class
	CellTissueTypeBasedGeneralisedLinearSpringForce<3>::OutputForceParameters(
			rParamsFile);
}

//// Explicit instantiation
//template class DirectionalAdhesionForce<1, 1> ;
//template class DirectionalAdhesionForce<1, 2> ;
//template class DirectionalAdhesionForce<2, 2> ;
//template class DirectionalAdhesionForce<1, 3> ;
//template class DirectionalAdhesionForce<2, 3> ;
//template class DirectionalAdhesionForce<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalAdhesionForce)
CHASTE_CLASS_EXPORT(DirectionalAdhesionForce)
