/*
 * OffLatticeSimulationDirectedDivision.hpp
 *
 *  Created on: May 16, 2017
 *      Author: Sonja Mathias
 */

#ifndef OFFLATTICESIMULATIONDIRECTEDDIVISION_HPP_
#define OFFLATTICESIMULATIONDIRECTEDDIVISION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "OffLatticeSimulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class OffLatticeSimulationDirectedDivision: public OffLatticeSimulation<
		ELEMENT_DIM, SPACE_DIM> {

	// Allow tests to access private members to test private functions
	friend class TestOffLatticeSimulationDirectedDivision;

private:

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive & boost::serialization::base_object<OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM> >(
						*this);
	}

	c_vector<double, SPACE_DIM> CalculateCellDivisionVector(CellPtr pParentCell);

public:

	OffLatticeSimulationDirectedDivision(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
			bool deleteCellPopulationInDestructor=false,
			bool initialiseCells=true);

	virtual ~OffLatticeSimulationDirectedDivision();

	void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulationDirectedDivision)

namespace boost {
namespace serialization {
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(Archive & ar,
		const OffLatticeSimulationDirectedDivision<ELEMENT_DIM, SPACE_DIM> * t,
		const unsigned int file_version) {
	// Save data required to construct instance
	const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population =
			&(t->rGetCellPopulation());
	ar & p_cell_population;
}

template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(Archive & ar,
		OffLatticeSimulationDirectedDivision<ELEMENT_DIM, SPACE_DIM> * t,
		const unsigned int file_version) {
	// Retrieve data from archive required to construct new instance
	AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_cell_population;
	ar & p_cell_population;

	// Invoke inplace constructor to initialise instance, last two variables set extra
	// member variables to be deleted as they are loaded from archive and to not initialise sells.
	::new (t) OffLatticeSimulationDirectedDivision<ELEMENT_DIM, SPACE_DIM>(*p_cell_population, true,
			false);
}
}
} // namespace

#endif /* OFFLATTICESIMULATIONDIRECTEDDIVISION_HPP_ */
