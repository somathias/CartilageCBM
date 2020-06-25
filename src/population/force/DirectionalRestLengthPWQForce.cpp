/*
 * DirectionalRestLengthPWQForce.cpp
 *
 *  Created on: June 24, 2020
 *      Author: Sonja Mathias
 */

#include "DirectionalRestLengthPWQForce.hpp"
#include "Debug.hpp"


DirectionalRestLengthPWQForce::DirectionalRestLengthPWQForce() :
		PWQGeneralisedLinearSpringForce<3>(){
}


c_vector<double, 3> DirectionalRestLengthPWQForce::CalculateForceBetweenNodes(
		unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
		AbstractCellPopulation<3>& rCellPopulation) {

	// We should only ever calculate the force between two distinct nodes
	assert(nodeAGlobalIndex != nodeBGlobalIndex);

	Node<3>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
	Node<3>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

	// Get the node locations
	c_vector<double, 3> node_a_location = p_node_a->rGetLocation();
	c_vector<double, 3> node_b_location = p_node_b->rGetLocation();

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

	/*
	 * If mUseCutOffLength has been set, then there is zero force between
	 * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
	 */
	if (this->mUseCutOffLength) {
		if (distance_between_nodes >= this->GetCutOffLength()) {
			return zero_vector<double>(3); // c_vector<double,3>() is not guaranteed to be fresh memory
		}
	}

	/*
	 * Calculate the rest length of the spring connecting the two nodes with a default
	 * value of 1.0.
	 */
	double rest_length_final = 1.0;

	if (bool(
			dynamic_cast<MeshBasedCellPopulation<3>*>(&rCellPopulation))) {
		rest_length_final = static_cast<MeshBasedCellPopulation<3,
				3>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex,
				nodeBGlobalIndex);
	} else if (bool(
			dynamic_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation))) {
		assert(node_a_radius > 0 && node_b_radius > 0);
		rest_length_final = node_a_radius + node_b_radius;
	}

	double rest_length = rest_length_final;

	CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(
			nodeAGlobalIndex);
	CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(
			nodeBGlobalIndex);

	double ageA = p_cell_A->GetAge();
	double ageB = p_cell_B->GetAge();

	assert(!std::isnan(ageA));
	assert(!std::isnan(ageB));

	/*
	 * If the cells are both newly divided, then the rest length of the spring
	 * connecting them grows linearly with time, until 1 hour after division.
	 */
	if (ageA
			< GeneralisedLinearSpringForce<3>::mMeinekeSpringGrowthDuration
			&& ageB
					< GeneralisedLinearSpringForce<3>::mMeinekeSpringGrowthDuration) {
		AbstractCentreBasedCellPopulation<3>* p_static_cast_cell_population =
				static_cast<AbstractCentreBasedCellPopulation<3,
						3>*>(&rCellPopulation);

		std::pair<CellPtr, CellPtr> cell_pair =
				p_static_cast_cell_population->CreateCellPair(p_cell_A,
						p_cell_B);

		if (p_static_cast_cell_population->IsMarkedSpring(cell_pair)) {
			// Spring rest length increases from a small value to the normal rest length over 1 hour
			double lambda =
					GeneralisedLinearSpringForce<3>::mMeinekeDivisionRestingSpringLength;
			rest_length = lambda
					+ (rest_length_final - lambda) * ageA
							/ GeneralisedLinearSpringForce<3,
									3>::mMeinekeSpringGrowthDuration;
		}
		if (ageA + SimulationTime::Instance()->GetTimeStep()
				>= GeneralisedLinearSpringForce<3>::mMeinekeSpringGrowthDuration) {
			// This spring is about to go out of scope
			p_static_cast_cell_population->UnmarkSpring(cell_pair);
		}
	}

	/*
	 * For apoptosis, progressively reduce the radius of the cell
	 */
	double a_rest_length = rest_length * 0.5;
	double b_rest_length = a_rest_length;

	if (bool(
			dynamic_cast<NodeBasedCellPopulation<3>*>(&rCellPopulation))) {
		assert(node_a_radius > 0 && node_b_radius > 0);
		a_rest_length = (node_a_radius / (node_a_radius + node_b_radius))
				* rest_length;
		b_rest_length = (node_b_radius / (node_a_radius + node_b_radius))
				* rest_length;
	}

	/*
	 * If either of the cells has begun apoptosis, then the length of the spring
	 * connecting them decreases linearly with time.
	 */
	if (p_cell_A->HasApoptosisBegun()) {
		double time_until_death_a = p_cell_A->GetTimeUntilDeath();
		a_rest_length = a_rest_length * time_until_death_a
				/ p_cell_A->GetApoptosisTime();
	}
	if (p_cell_B->HasApoptosisBegun()) {
		double time_until_death_b = p_cell_B->GetTimeUntilDeath();
		b_rest_length = b_rest_length * time_until_death_b
				/ p_cell_B->GetApoptosisTime();
	}

	rest_length = a_rest_length + b_rest_length;

    /*
     * Scale with the direction of the cell-cell interaction. Calculate absolute value of
     * scalar product with normal to x-y plane.
     * 
     */
    double abs_dot_product  = fabs(unit_difference[2]);
    if(abs_dot_product < 1.0/sqrt(2.0)){
        rest_length = ((2.0-sqrt(2.0))*abs_dot_product+1.0)*rest_length;
    }
    else{
        rest_length = (1.0/(1-sqrt(2.0))*((2.0-sqrt(2.0))*abs_dot_product-1.0))*rest_length;
    }

	//assert(rest_length <= 1.0+1e-12); ///\todo #1884 Magic number: would "<= 1.0" do?

	// Although in this class the 'spring constant' is a constant parameter, in
	// subclasses it can depend on properties of each of the cells
	double overlap = distance_between_nodes - rest_length;
	bool is_closer_than_rest_length = (overlap <= 0);
	double multiplication_factor = this->VariableSpringConstantMultiplicationFactor(
			nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation,
			is_closer_than_rest_length);
	double spring_stiffness_adhesion = GeneralisedLinearSpringForce<3,
			3>::mMeinekeSpringStiffness;
    double spring_stiffness_repulsion = mRepulsionSpringStiffness;

    double cut_off_length = AbstractTwoBodyInteractionForce<3>::mMechanicsCutOffLength;
    double ratio = spring_stiffness_adhesion / spring_stiffness_repulsion;
    double cut_off_repulsion = rest_length / (1 - sqrt(ratio) * (1-rest_length / cut_off_length));

	if (bool(
			dynamic_cast<MeshBasedCellPopulation<3>*>(&rCellPopulation))) {
		return multiplication_factor * spring_stiffness_adhesion
				* unit_difference * overlap;
	} else {
        if (distance_between_nodes < cut_off_repulsion)
		{
			c_vector<double, 3> temp = multiplication_factor
				    * unit_difference 
                    * (spring_stiffness_adhesion 
                        * (1 - distance_between_nodes / cut_off_length) 
                        * (1 - distance_between_nodes / cut_off_length)
                       - spring_stiffness_repulsion 
                        * (1 - distance_between_nodes / cut_off_repulsion) 
                        * (1 - distance_between_nodes / cut_off_repulsion)
                      );
			return temp;
		} else {
			c_vector<double, 3> temp = multiplication_factor
				    * unit_difference 
                    * spring_stiffness_adhesion 
                    * (1 - distance_between_nodes / cut_off_length) 
                    * (1 - distance_between_nodes / cut_off_length);
			return temp;
		}

	}
}


void DirectionalRestLengthPWQForce::OutputForceParameters(out_stream& rParamsFile) {

	// Call direct parent class
	PWQGeneralisedLinearSpringForce<3>::OutputForceParameters(
			rParamsFile);
}

// // Explicit instantiation
// template class DirectionalRestLengthPWQForce<1, 1> ;
// template class DirectionalRestLengthPWQForce<1, 2> ;
// template class DirectionalRestLengthPWQForce<2, 2> ;
// template class DirectionalRestLengthPWQForce<1, 3> ;
// template class DirectionalRestLengthPWQForce<2, 3> ;
// template class DirectionalRestLengthPWQForce<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(DirectionalRestLengthPWQForce)
CHASTE_CLASS_EXPORT(DirectionalRestLengthPWQForce)

