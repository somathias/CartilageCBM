#ifndef REPULSIONCUBICFORCE_HPP_
#define REPULSIONCUBICFORCE_HPP_

#include "CubicGeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"

template<unsigned DIM>
class RepulsionCubicForce : public CubicGeneralisedLinearSpringForce<DIM>
{
private :

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<CubicGeneralisedLinearSpringForce<DIM> >(*this);
    }

public :

    RepulsionCubicForce();

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionCubicForce)

#endif /*REPULSIONCUBICFORCE_HPP_*/
