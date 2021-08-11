#include "OrientationBasedDivisionRule.hpp"
#include "RandomNumberGenerator.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "HorizontalCellDivisionDirection.hpp"
#include "FixedCellDivisionDirection.hpp"
#include "PerichondrialCellTissueType.hpp"


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM>> OrientationBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation)
{
    // Get separation parameter
    double separation = rCellPopulation.GetMeinekeDivisionSeparation();

    if (pParentCell->HasCellProperty<PerichondrialCellTissueType>()){
        //parent cell is perichondrial - do not place cells at overlap
        separation = 1.0;
    }

    // Make a random direction vector of the required length
    c_vector<double, SPACE_DIM> main_direction;

    /*
	 * Check if the cell has a the property AbstractCellDivisionDirection
	 * We need to check both cases independently.
	 */
    if (pParentCell->HasCellProperty<UpwardsCellDivisionDirection<SPACE_DIM>>())
    {
        /*
		 * Use the main cell division direction provided
		 */
        CellPropertyCollection collection =
            pParentCell->rGetCellPropertyCollection();
        CellPropertyCollection direction_collection =
            collection.GetProperties<
                UpwardsCellDivisionDirection<SPACE_DIM>>();

        if (direction_collection.GetSize() == 1)
        {
            boost::shared_ptr<UpwardsCellDivisionDirection<SPACE_DIM>> p_direction =
                boost::static_pointer_cast<
                    UpwardsCellDivisionDirection<SPACE_DIM>>(
                    direction_collection.GetProperty());
            main_direction = 0.5 * separation * p_direction->GetCellDivisionDirection();
        }
    }
    else if (pParentCell->HasCellProperty<
                 DownwardsCellDivisionDirection<SPACE_DIM>>())
    {
        /*
		 * Use the main cell division direction provided
		 */
        CellPropertyCollection collection =
            pParentCell->rGetCellPropertyCollection();
        CellPropertyCollection direction_collection =
            collection.GetProperties<
                DownwardsCellDivisionDirection<SPACE_DIM>>();

        if (direction_collection.GetSize() == 1)
        {
            boost::shared_ptr<DownwardsCellDivisionDirection<SPACE_DIM>> p_direction =
                boost::static_pointer_cast<
                    DownwardsCellDivisionDirection<SPACE_DIM>>(
                    direction_collection.GetProperty());
            main_direction = 0.5 * separation * p_direction->GetCellDivisionDirection();
        }
    }
    else if(pParentCell->HasCellProperty<
                 HorizontalCellDivisionDirection<SPACE_DIM>>())
    {

        CellPropertyCollection collection =
            pParentCell->rGetCellPropertyCollection();
        CellPropertyCollection direction_collection =
            collection.GetProperties<
                HorizontalCellDivisionDirection<SPACE_DIM>>();

        if (direction_collection.GetSize() == 1)
        {
            double random_angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

            main_direction(0) = 0.5 * separation * cos(random_angle);
            main_direction(1) = 0.5 * separation * sin(random_angle);
        }
    }
    else if(pParentCell->HasCellProperty<
                 FixedCellDivisionDirection<SPACE_DIM>>())
    {

        CellPropertyCollection collection =
            pParentCell->rGetCellPropertyCollection();
        CellPropertyCollection direction_collection =
            collection.GetProperties<
                FixedCellDivisionDirection<SPACE_DIM>>();

        if (direction_collection.GetSize() == 1)
        {
            boost::shared_ptr<FixedCellDivisionDirection<SPACE_DIM>> p_direction =
                boost::static_pointer_cast<
                    FixedCellDivisionDirection<SPACE_DIM>>(
                    direction_collection.GetProperty());
            main_direction = 0.5 * separation * p_direction->GetCellDivisionDirection();
        }
    }
    else
    {
        /*
		 * Pick a random direction and move the parent cell backwards by 0.5*separation
		 * in that direction and return the position of the daughter cell 0.5*separation
		 * forwards in that direction.
		 */
        switch (SPACE_DIM)
        {
        case 1:
        {
            double random_direction = -1.0 + 2.0 * (RandomNumberGenerator::Instance()->ranf() < 0.5);

            main_direction(0) = 0.5 * separation * random_direction;
            break;
        }
        case 2:
        {
            double random_angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();

            main_direction(0) = 0.5 * separation * cos(random_angle);
            main_direction(1) = 0.5 * separation * sin(random_angle);
            break;
        }
        case 3:
        {
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

            main_direction(0) = 0.5 * separation * cos(random_azimuth_angle) * sin(random_zenith_angle);
            main_direction(1) = 0.5 * separation * sin(random_azimuth_angle) * sin(random_zenith_angle);
            main_direction(2) = 0.5 * separation * cos(random_zenith_angle);
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
        }
    }

    c_vector<double, SPACE_DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    if (!pParentCell->HasCellProperty<PerichondrialCellTissueType>()){
        //parent cell is not perichondrial - move backwards
        parent_position -=  main_direction;
    }

    c_vector<double, SPACE_DIM> daughter_position = parent_position + 2 * main_direction;

    std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM>> positions(parent_position, daughter_position);

    return positions;
}

// Explicit instantiation
template class OrientationBasedDivisionRule<1, 1>;
template class OrientationBasedDivisionRule<1, 2>;
template class OrientationBasedDivisionRule<2, 2>;
template class OrientationBasedDivisionRule<1, 3>;
template class OrientationBasedDivisionRule<2, 3>;
template class OrientationBasedDivisionRule<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OrientationBasedDivisionRule)
