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
				1.0), mAlpha(5.0) {
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

/**
 * Set mAlpha.
 *
 * @param alpha the new value of mAlpha
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::SetAlpha(
		double alpha) {
	assert(alpha > 0.0);
	mAlpha = alpha;
}

/**
 * @return #mAlpha
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypeBasedGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::GetAlpha() {
	return mAlpha;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> CellTissueTypeBasedGeneralisedLinearSpringForce<
		ELEMENT_DIM, SPACE_DIM>::CalculateForceBetweenNodes(
		unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation) {
	// We should only ever calculate the force between two distinct nodes
	assert(nodeAGlobalIndex != nodeBGlobalIndex);

	Node<SPACE_DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
	Node<SPACE_DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

	// Get the node locations
	c_vector<double, SPACE_DIM> node_a_location = p_node_a->rGetLocation();
	c_vector<double, SPACE_DIM> node_b_location = p_node_b->rGetLocation();

	// Get the node radii for a NodeBasedCellPopulation
	double node_a_radius = 0.0;
	double node_b_radius = 0.0;

	if (bool(
			dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation))) {
		node_a_radius = p_node_a->GetRadius();
		node_b_radius = p_node_b->GetRadius();
	}

	// Get the unit vector parallel to the line joining the two nodes
	c_vector<double, SPACE_DIM> unit_difference;
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
			return zero_vector<double>(SPACE_DIM); // c_vector<double,SPACE_DIM>() is not guaranteed to be fresh memory
		}
	}

	/*
	 * Calculate the rest length of the spring connecting the two nodes with a default
	 * value of 1.0.
	 */
	double rest_length_final = 1.0;

	if (bool(
			dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation))) {
		rest_length_final = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,
				SPACE_DIM>*>(&rCellPopulation)->GetRestLength(nodeAGlobalIndex,
				nodeBGlobalIndex);
	} else if (bool(
			dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation))) {
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
			< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::mMeinekeSpringGrowthDuration
			&& ageB
					< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::mMeinekeSpringGrowthDuration) {
		AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_static_cast_cell_population =
				static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,
						SPACE_DIM>*>(&rCellPopulation);

		std::pair<CellPtr, CellPtr> cell_pair =
				p_static_cast_cell_population->CreateCellPair(p_cell_A,
						p_cell_B);

		if (p_static_cast_cell_population->IsMarkedSpring(cell_pair)) {
			// Spring rest length increases from a small value to the normal rest length over 1 hour
			double lambda =
					GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::mMeinekeDivisionRestingSpringLength;
			rest_length = lambda
					+ (rest_length_final - lambda) * ageA
							/ GeneralisedLinearSpringForce<ELEMENT_DIM,
									SPACE_DIM>::mMeinekeSpringGrowthDuration;
		}
		if (ageA + SimulationTime::Instance()->GetTimeStep()
				>= GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::mMeinekeSpringGrowthDuration) {
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
			dynamic_cast<NodeBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation))) {
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
	//assert(rest_length <= 1.0+1e-12); ///\todo #1884 Magic number: would "<= 1.0" do?

	// Although in this class the 'spring constant' is a constant parameter, in
	// subclasses it can depend on properties of each of the cells
	double overlap = distance_between_nodes - rest_length;
	bool is_closer_than_rest_length = (overlap <= 0);
	double multiplication_factor = VariableSpringConstantMultiplicationFactor(
			nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation,
			is_closer_than_rest_length);
	double spring_stiffness = GeneralisedLinearSpringForce<ELEMENT_DIM,
			SPACE_DIM>::mMeinekeSpringStiffness;

	if (bool(
			dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(&rCellPopulation))) {
		return multiplication_factor * spring_stiffness * unit_difference
				* overlap;
	} else {
		// A reasonably stable simple force law
		if (is_closer_than_rest_length) //overlap is negative
		{
			//log(x+1) is undefined for x<=-1
			assert(overlap > -rest_length_final);
			c_vector<double, SPACE_DIM> temp = multiplication_factor
					* spring_stiffness * unit_difference * rest_length_final
					* log(1.0 + overlap / rest_length_final);
			return temp;
		} else {
			//double alpha = 5.0; we want to use our member variable instead
			c_vector<double, SPACE_DIM> temp = multiplication_factor
					* spring_stiffness * unit_difference * overlap
					* exp(-mAlpha * overlap / rest_length_final);
			return temp;
		}
	}
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

