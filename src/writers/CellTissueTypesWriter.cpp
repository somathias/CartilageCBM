#include "CellTissueTypesWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellTissueTypesWriter<ELEMENT_DIM, SPACE_DIM>::CellTissueTypesWriter() :
		AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizcelltissuetypes") {
	this->mVtkCellDataName = "Cell tissue types";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTissueTypesWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
		CellPtr pCell,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) {
	double colour = 0.0;
	if (pCell->HasCellProperty<PerichondrialCellTissueType>()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				PerichondrialCellTissueType>();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<PerichondrialCellTissueType> p_tissue_type =
				boost::static_pointer_cast<PerichondrialCellTissueType>(
						direction_collection.GetProperty());

		colour = p_tissue_type->GetColour();

	} else if (pCell->HasCellProperty<ChondrocyteCellTissueType>()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				ChondrocyteCellTissueType>();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<ChondrocyteCellTissueType> p_tissue_type =
				boost::static_pointer_cast<ChondrocyteCellTissueType>(
						direction_collection.GetProperty());

		colour = p_tissue_type->GetColour();
	}
	return colour;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTissueTypesWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell,
		AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) {
	unsigned colour = 0;
	if (pCell->HasCellProperty<PerichondrialCellTissueType>()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				PerichondrialCellTissueType>();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<PerichondrialCellTissueType> p_tissue_type =
				boost::static_pointer_cast<PerichondrialCellTissueType>(
						direction_collection.GetProperty());

		colour = p_tissue_type->GetColour();

	} else if (pCell->HasCellProperty<ChondrocyteCellTissueType>()) {

		CellPropertyCollection collection = pCell->rGetCellPropertyCollection();
		CellPropertyCollection direction_collection = collection.GetProperties<
				ChondrocyteCellTissueType>();

		assert(direction_collection.GetSize() == 1);
		boost::shared_ptr<ChondrocyteCellTissueType> p_tissue_type =
				boost::static_pointer_cast<ChondrocyteCellTissueType>(
						direction_collection.GetProperty());

		colour = p_tissue_type->GetColour();
	}

	*this->mpOutStream << colour << " ";
}

// Explicit instantiation
template class CellTissueTypesWriter<1, 1> ;
template class CellTissueTypesWriter<1, 2> ;
template class CellTissueTypesWriter<2, 2> ;
template class CellTissueTypesWriter<1, 3> ;
template class CellTissueTypesWriter<2, 3> ;
template class CellTissueTypesWriter<3, 3> ;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellTissueTypesWriter)
