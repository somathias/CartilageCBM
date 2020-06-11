/*
 * DirectionalRepulsionPWQForce.cpp
 *
 *  Created on: June 10, 2020
 *      Author: Sonja Mathias
 */

#include "DirectionalRepulsionPWQForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"
#include <iostream>

DirectionalRepulsionPWQForce::DirectionalRepulsionPWQForce() :
		PWQGeneralisedLinearSpringForce<3>(), mBaselineRepulsionMultiplier(
				0.1) {

}

double DirectionalRepulsionPWQForce::GetBaselineRepulsionMultiplier() {
	return mBaselineRepulsionMultiplier;
}


void DirectionalRepulsionPWQForce::SetBaselineRepulsionMultiplier(
		double baselineRepulsionMultiplier) {
	assert(baselineRepulsionMultiplier > 0.0);
	mBaselineRepulsionMultiplier =
			baselineRepulsionMultiplier;
}

double DirectionalRepulsionPWQForce::VariableSpringConstantMultiplicationFactor(
		unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
		AbstractCellPopulation<3>& rCellPopulation,
		bool isCloserThanRestLength) {

	if (!isCloserThanRestLength) {
		return 1.0;
	} else {

		Node<3>* p_node_a = rCellPopulation.GetNode(
					nodeAGlobalIndex);
		Node<3>* p_node_b = rCellPopulation.GetNode(
					nodeBGlobalIndex);

		// Get the node locations
		c_vector<double, 3> node_a_location =
					p_node_a->rGetLocation();
		c_vector<double, 3> node_b_location =
					p_node_b->rGetLocation();


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


		double abs_dot_product  = fabs(unit_difference[2]);
        double directional_multiplier = abs_dot_product * (1.0-mBaselineRepulsionMultiplier) + mBaselineRepulsionMultiplier;
        PRINT_VARIABLE(abs_dot_product);
        PRINT_VARIABLE(directional_multiplier);
        return directional_multiplier;
	}
}

void DirectionalRepulsionPWQForce::OutputForceParameters(out_stream& rParamsFile) {

	*rParamsFile << "\t\t\t<BaselineRepulsionMultiplier>"
			<< mBaselineRepulsionMultiplier << "</BaselineRepulsionMultiplier>\n";

	// Call direct parent class
	PWQGeneralisedLinearSpringForce<3>::OutputForceParameters(
			rParamsFile);
}

//// Explicit instantiation
//template class DirectionalRepulsionPWQForce<1, 1> ;
//template class DirectionalRepulsionPWQForce<1, 2> ;
//template class DirectionalRepulsionPWQForce<2, 2> ;
//template class DirectionalRepulsionPWQForce<1, 3> ;
//template class DirectionalRepulsionPWQForce<2, 3> ;
//template class DirectionalRepulsionPWQForce<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalRepulsionPWQForce)
CHASTE_CLASS_EXPORT(DirectionalRepulsionPWQForce)
