#ifndef OFFLATTICESIMULATION2DDIRECTEDDIVISION_HPP_
#define OFFLATTICESIMULATION2DDIRECTEDDIVISION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "OffLatticeSimulation.hpp"

/* I based the structure of this file on CryptSimulation1d */

class OffLatticeSimulation2dDirectedDivision : public OffLatticeSimulation<2>
{
    // Allow tests to access private members to test private functions
    friend class TestOffLatticeSimulation2dDirectedDivision;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<2> >(*this);
    }

    c_vector<double, 2> CalculateCellDivisionVector(CellPtr pParentCell);

public:

    OffLatticeSimulation2dDirectedDivision(AbstractCellPopulation<2>& rCellPopulation,
                      bool deleteCellPopulationInDestructor=false,
                      bool initialiseCells=true);

    virtual ~OffLatticeSimulation2dDirectedDivision();

    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulation2dDirectedDivision)

namespace boost
{
namespace serialization
{
template<class Archive>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulation2dDirectedDivision * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

template<class Archive>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulation2dDirectedDivision * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar & p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)OffLatticeSimulation2dDirectedDivision(*p_cell_population, true, false);
}
}
} // namespace

#endif /*OFFLATTICESIMULATION2DDIRECTEDDIVISION_HPP_*/
