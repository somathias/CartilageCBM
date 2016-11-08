#include "OffLatticeSimulation2dDirectedDivision.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"

/* 
 * Implements a directed cell division.
 *
 * I based the structure of this file on CryptSimulation1d */

OffLatticeSimulation2dDirectedDivision::OffLatticeSimulation2dDirectedDivision(
  AbstractCellPopulation<2>& rCellPopulation, 
  bool deleteCellPopulationInDestructor, 
  bool initialiseCells) 
  : OffLatticeSimulation<2>(rCellPopulation,
			    deleteCellPopulationInDestructor,
			    initialiseCells)
{
}

OffLatticeSimulation2dDirectedDivision::~OffLatticeSimulation2dDirectedDivision()
{
}

c_vector<double, 2> OffLatticeSimulation2dDirectedDivision::CalculateCellDivisionVector(CellPtr pParentCell)
{
     // If it is not vertex based...
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<2>*>(&(this->mrCellPopulation))))
    {
        // Location of parent and daughter cells
        c_vector<double, 2> parent_coords = this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, 2> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<2>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

        // Make a main direction vector of the required length
        c_vector<double, 2> main_direction = zero_vector<double>(2);
        main_direction(1) = 1.0;          //upwards
	
        c_vector<double, 2> new_parent_coords = parent_coords - 0.5*separation*main_direction;
        daughter_coords = parent_coords + 0.5*separation*main_direction;
	
        // Set the parent to use this location
        ChastePoint<2> parent_coords_point(new_parent_coords);
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        this->mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else if (bool(dynamic_cast<VertexBasedCellPopulation<2>*>(&(this->mrCellPopulation))))
    {
        VertexBasedCellPopulation<2>* p_vertex_population = dynamic_cast<VertexBasedCellPopulation<2>*>(&(this->mrCellPopulation));
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule = p_vertex_population->GetVertexBasedDivisionRule();

        return p_division_rule->CalculateCellDivisionVector(pParentCell, *p_vertex_population);
    }
    else
    {
        // All classes derived from AbstractOffLatticeCellPopulation are covered by the above (except user-derived classes),
        // i.e. if you want to use this class with your own subclass of AbstractOffLatticeCellPopulation, then simply
        // comment out the line below
        NEVER_REACHED;
    }
}

void OffLatticeSimulation2dDirectedDivision::OutputSimulationParameters(out_stream& rParamsFile)
{
    // No parameters to output

    // Call method on direct parent class
    OffLatticeSimulation<2>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulation2dDirectedDivision)
