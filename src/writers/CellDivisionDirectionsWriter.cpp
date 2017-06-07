#include "CellDivisionDirectionsWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellDivisionDirectionsWriter<ELEMENT_DIM, SPACE_DIM>::CellDivisionDirectionsWriter() :
		AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>(
				"results.vizcelldivisiondirections") {
	this->mVtkCellDataName = "Cell division directions";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellDivisionDirectionsWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
		CellPtr pCell,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) {

	double colour = 0.0;
	if (pCell->HasCellProperty<UpwardsCellDivisionDirection<SPACE_DIM> >()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				UpwardsCellDivisionDirection<SPACE_DIM> >();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<UpwardsCellDivisionDirection<SPACE_DIM> > p_direction =
				boost::static_pointer_cast<
						UpwardsCellDivisionDirection<SPACE_DIM> >(
						direction_collection.GetProperty());

		colour = p_direction->GetColour();

	} else if (pCell->HasCellProperty<DownwardsCellDivisionDirection<SPACE_DIM> >()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				DownwardsCellDivisionDirection<SPACE_DIM> >();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<DownwardsCellDivisionDirection<SPACE_DIM> > p_direction =
				boost::static_pointer_cast<
						DownwardsCellDivisionDirection<SPACE_DIM> >(
						direction_collection.GetProperty());

		colour = p_direction->GetColour();
	}
	return colour;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellDivisionDirectionsWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(
		CellPtr pCell,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) {

	unsigned colour = 0;

	if (pCell->HasCellProperty<UpwardsCellDivisionDirection<SPACE_DIM> >()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				UpwardsCellDivisionDirection<SPACE_DIM> >();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<UpwardsCellDivisionDirection<SPACE_DIM> > p_direction =
				boost::static_pointer_cast<
						UpwardsCellDivisionDirection<SPACE_DIM> >(
						direction_collection.GetProperty());

		colour = p_direction->GetColour();

	} else if (pCell->HasCellProperty<DownwardsCellDivisionDirection<SPACE_DIM> >()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				DownwardsCellDivisionDirection<SPACE_DIM> >();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<DownwardsCellDivisionDirection<SPACE_DIM> > p_direction =
				boost::static_pointer_cast<
						DownwardsCellDivisionDirection<SPACE_DIM> >(
						direction_collection.GetProperty());

		colour = p_direction->GetColour();
	}

	*this->mpOutStream << colour << " ";
}

// Explicit instantiation
template class CellDivisionDirectionsWriter<1, 1> ;
template class CellDivisionDirectionsWriter<1, 2> ;
template class CellDivisionDirectionsWriter<2, 2> ;
template class CellDivisionDirectionsWriter<1, 3> ;
template class CellDivisionDirectionsWriter<2, 3> ;
template class CellDivisionDirectionsWriter<3, 3> ;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellDivisionDirectionsWriter)
