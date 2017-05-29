/*
 * OffLatticeSimulationDirectedDivision.cpp
 *
 *  Created on: May 16, 2017
 *      Author: Sonja Mathias
 */
#include "OffLatticeSimulationDirectedDivision.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"

#include "AbstractCellDivisionDirection.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"

/*
 * Implements a directed cell division.
 *
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> OffLatticeSimulationDirectedDivision<
		ELEMENT_DIM, SPACE_DIM>::OffLatticeSimulationDirectedDivision(
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
		bool deleteCellPopulationInDestructor, bool initialiseCells) :
		OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>(rCellPopulation,
				deleteCellPopulationInDestructor, initialiseCells) {
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> OffLatticeSimulationDirectedDivision<
		ELEMENT_DIM, SPACE_DIM>::~OffLatticeSimulationDirectedDivision() {
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> c_vector<double, SPACE_DIM> OffLatticeSimulationDirectedDivision<
		ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
		CellPtr pParentCell) {
	// If it is not vertex based...
	if (bool(
			dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,
					SPACE_DIM>*>(&(this->mrCellPopulation)))) {
		// Location of parent and daughter cells
		c_vector<double, SPACE_DIM> parent_coords =
				this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
		c_vector<double, SPACE_DIM> daughter_coords;

		// Get separation parameter
		double separation =
				static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,
						SPACE_DIM>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

		c_vector<double, SPACE_DIM> main_direction;

		/*
		 * Check if the cell has a the property AbstractCellDivisionDirection
		 * We need to check both cases independently.
		 */
		if (pParentCell->HasCellProperty<UpwardsCellDivisionDirection<SPACE_DIM> >()) {
			/*
			 * Use the main cell division direction provided
			 */
			CellPropertyCollection collection =
					pParentCell->rGetCellPropertyCollection();
			CellPropertyCollection direction_collection =
					collection.GetProperties<
							UpwardsCellDivisionDirection<SPACE_DIM> >();

			if (direction_collection.GetSize() == 1) {
				boost::shared_ptr<UpwardsCellDivisionDirection<SPACE_DIM> > p_direction =
						boost::static_pointer_cast<
								UpwardsCellDivisionDirection<SPACE_DIM> >(
								direction_collection.GetProperty());
				main_direction = 0.5 * separation
						* p_direction->GetCellDivisionDirection();
			}
		} else if (pParentCell->HasCellProperty<
				DownwardsCellDivisionDirection<SPACE_DIM> >()) {
			/*
			 * Use the main cell division direction provided
			 */
			CellPropertyCollection collection =
					pParentCell->rGetCellPropertyCollection();
			CellPropertyCollection direction_collection =
					collection.GetProperties<
							DownwardsCellDivisionDirection<SPACE_DIM> >();

			if (direction_collection.GetSize() == 1) {
				boost::shared_ptr<DownwardsCellDivisionDirection<SPACE_DIM> > p_direction =
						boost::static_pointer_cast<
								DownwardsCellDivisionDirection<SPACE_DIM> >(
								direction_collection.GetProperty());
				main_direction = 0.5 * separation
						* p_direction->GetCellDivisionDirection();
			}
		} else {
			/*
			 * Pick a random direction and move the parent cell backwards by 0.5*separation
			 * in that direction and return the position of the daughter cell 0.5*separation
			 * forwards in that direction.
			 */
			switch (SPACE_DIM) {
			case 1: {
				double random_direction = -1.0
						+ 2.0
								* (RandomNumberGenerator::Instance()->ranf()
										< 0.5);

				main_direction(0) = 0.5 * separation * random_direction;
				break;
			}
			case 2: {
				double random_angle = 2.0 * M_PI
						* RandomNumberGenerator::Instance()->ranf();

				main_direction(0) = 0.5 * separation * cos(random_angle);
				main_direction(1) = 0.5 * separation * sin(random_angle);
				break;
			}
			case 3: {
				/*
				 * Note that to pick a random point on the surface of a sphere, it is incorrect
				 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
				 * [0, pi) respectively, since points picked in this way will be 'bunched' near
				 * the poles. See #2230.
				 */
				double u = RandomNumberGenerator::Instance()->ranf();
				double v = RandomNumberGenerator::Instance()->ranf();

				double random_azimuth_angle = 2 * M_PI * u;
				double random_zenith_angle = std::acos(2 * v - 1);

				main_direction(0) = 0.5 * separation * cos(random_azimuth_angle)
						* sin(random_zenith_angle);
				main_direction(1) = 0.5 * separation * sin(random_azimuth_angle)
						* sin(random_zenith_angle);
				main_direction(2) = 0.5 * separation * cos(random_zenith_angle);
				break;
			}
			default:
				// This can't happen
				NEVER_REACHED;
				}
			}
//			// Make a main direction vector of the required length
//		c_vector<double, SPACE_DIM> main_direction = zero_vector<double>(
//				SPACE_DIM);
//		main_direction(SPACE_DIM - 1) = 1.0;          //upwards

		c_vector<double, SPACE_DIM> new_parent_coords = parent_coords
				- main_direction;
		daughter_coords = parent_coords + main_direction;

		// Set the parent to use this location
		ChastePoint<SPACE_DIM> parent_coords_point(new_parent_coords);
		unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(
				pParentCell);
		this->mrCellPopulation.SetNode(node_index, parent_coords_point);

		return daughter_coords;
	} else if (bool(
			dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation)))) {
		VertexBasedCellPopulation<SPACE_DIM> *p_vertex_population =
				dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation));
		boost::shared_ptr<AbstractVertexBasedDivisionRule<SPACE_DIM> > p_division_rule =
				p_vertex_population->GetVertexBasedDivisionRule();

		return p_division_rule->CalculateCellDivisionVector(pParentCell,
				*p_vertex_population);
	} else {
		// All classes derived from AbstractOffLatticeCellPopulation are covered by the above (except user-derived classes),
		// i.e. if you want to use this class with your own subclass of AbstractOffLatticeCellPopulation, then simply
		// comment out the line below
		NEVER_REACHED;
		}
	}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> void OffLatticeSimulationDirectedDivision<
		ELEMENT_DIM, SPACE_DIM>::OutputSimulationParameters(
		out_stream& rParamsFile) {
	// No parameters to output

	// Call method on direct parent class
	OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>::OutputSimulationParameters(
			rParamsFile);
}

template class OffLatticeSimulationDirectedDivision<1, 1> ;
template class OffLatticeSimulationDirectedDivision<1, 2> ;
template class OffLatticeSimulationDirectedDivision<2, 2> ;
template class OffLatticeSimulationDirectedDivision<1, 3> ;
template class OffLatticeSimulationDirectedDivision<2, 3> ;
template class OffLatticeSimulationDirectedDivision<3, 3> ;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulationDirectedDivision)
