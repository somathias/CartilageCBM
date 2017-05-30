/*
 * CellTissueTypeBasedCellCycleModel.cpp
 *
 *  Created on: May 29, 2017
 *      Author: Sonja Mathias
 */

#include "CellTissueTypeBasedCellCycleModel.hpp"

CellTissueTypeBasedCellCycleModel::CellTissueTypeBasedCellCycleModel() {

}

AbstractCellCycleModel* CellTissueTypeBasedCellCycleModel::CreateCellCycleModel() {
	// Create a new cell-cycle model
	CellTissueTypeBasedCellCycleModel* p_model =
			new CellTissueTypeBasedCellCycleModel();

	/*
	 * Set each member variable of the new cell-cycle model that inherits
	 * its value from the parent.
	 *
	 * Note 1: some of the new cell-cycle model's member variables (namely
	 * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
	 * correctly initialized in its constructor.
	 *
	 * Note 2: one or more of the new cell-cycle model's member variables
	 * may be set/overwritten as soon as InitialiseDaughterCell() is called on
	 * the new cell-cycle model.
	 *
	 * Note 3: the member variable mDimension remains unset, since this cell-cycle
	 * model does not need to know the spatial dimension, so if we were to call
	 * SetDimension() on the new cell-cycle model an exception would be triggered;
	 * hence we do not set this member variable.
	 */
	p_model->SetBirthTime(mBirthTime);
	p_model->SetMinimumGapDuration(mMinimumGapDuration);
	p_model->SetStemCellG1Duration(mStemCellG1Duration);
	p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
	p_model->SetSDuration(mSDuration);
	p_model->SetG2Duration(mG2Duration);
	p_model->SetMDuration(mMDuration);
	p_model->SetGeneration(mGeneration);
	p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

	return p_model;
}

void CellTissueTypeBasedCellCycleModel::InitialiseDaughterCell() {

	/*
	 * In this cell tissue type based cell-cycle model, the daughter cell
	 * of a perichondrial cell is of type chondrocyte.
	 */

	if (mpCell->HasCellProperty<PerichondrialCellTissueType>()) {

		//Get old cell tissue type and delete it
		CellPropertyCollection tissue_type_collection =
				mpCell->rGetCellPropertyCollection().GetProperties<
				PerichondrialCellTissueType>();

		assert(tissue_type_collection.GetSize() == 1);

		boost::shared_ptr<AbstractCellProperty> p_old_type = tissue_type_collection.GetProperty();
//		boost::shared_ptr<PerichondrialCellTissueType> p_old_tissue_type =
//				boost::static_pointer_cast<PerichondrialCellTissueType>(p_old_type);

		p_old_type->DecrementCellCount();
		mpCell->rGetCellPropertyCollection().RemoveProperty(p_old_type);

		//Add new chondrocyte type
		boost::shared_ptr<AbstractCellProperty> p_chondrocyte_type =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>();

		mpCell->AddCellProperty(p_chondrocyte_type);
	}

	StochasticDurationGenerationBasedCellCycleModel::InitialiseDaughterCell();
}

void CellTissueTypeBasedCellCycleModel::OutputCellCycleModelParameters(
		out_stream& rParamsFile) {
	// No new parameters to output

	// Call method on direct parent class
	StochasticDurationGenerationBasedCellCycleModel::OutputCellCycleModelParameters(
			rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellTissueTypeBasedCellCycleModel)
