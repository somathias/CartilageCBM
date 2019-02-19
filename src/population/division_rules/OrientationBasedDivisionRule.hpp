#ifndef ORIENTATIONBASEDDIVISIONRULE_HPP_
#define ORIENTATIONBASEDDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

#include "UblasVectorInclude.hpp" //needed for c_vector reference

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class OrientationBasedDivisionRule : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    OrientationBasedDivisionRule()
    {
    }

    virtual ~OrientationBasedDivisionRule()
    {
    }

    virtual std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OrientationBasedDivisionRule)

#endif // ORIENTATIONBASEDDIVISIONRULE_HPP_