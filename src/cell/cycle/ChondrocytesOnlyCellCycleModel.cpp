/*
 * ChondrocytesOnlyCellCycleModel.cpp
 *
 *  Created on: May 29, 2017
 *      Author: Sonja Mathias
 */

#include "ChondrocytesOnlyCellCycleModel.hpp"
#include <iostream>


ChondrocytesOnlyCellCycleModel::ChondrocytesOnlyCellCycleModel() 
	: AbstractSimpleGenerationalCellCycleModel(), mPatchSizeLimit(6) {

}

AbstractCellCycleModel* ChondrocytesOnlyCellCycleModel::CreateCellCycleModel() {
	// Create a new cell-cycle model
	ChondrocytesOnlyCellCycleModel* p_model =
			new ChondrocytesOnlyCellCycleModel();

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

	p_model->SetPatchSizeLimit(mPatchSizeLimit);

	return p_model;
}

void ChondrocytesOnlyCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    double uniform_random_number = RandomNumberGenerator::Instance()->ranf();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = -log(uniform_random_number)*GetStemCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = -log(uniform_random_number)*GetTransitCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}


bool ChondrocytesOnlyCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ((mCurrentCellCyclePhase != G_ZERO_PHASE) &&
            (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()))
        {
            mReadyToDivide = true;
        }
    }
	//check size of clonal patch before dividing - this fails if "patch size" has not been set - not sure why
	try {
		if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() && mpCell->GetCellData()->GetItem("patch size") >= mPatchSizeLimit){
		//if (mpCell->GetCellData()->GetItem("patch size") >= mPatchSizeLimit){
			mReadyToDivide = false;
		}
	} catch (const std::exception& e){
		std::cout << "PatchSizeTrackingModifier has not been enabled." << std::endl;
	}
    return mReadyToDivide;
}

void ChondrocytesOnlyCellCycleModel::SetPatchSizeLimit(unsigned patchSizeLimit){
	assert(patchSizeLimit >= 1);

	mPatchSizeLimit = patchSizeLimit;
}

unsigned ChondrocytesOnlyCellCycleModel::GetPatchSizeLimit(){
	return mPatchSizeLimit;
}


void ChondrocytesOnlyCellCycleModel::OutputCellCycleModelParameters(
		out_stream& rParamsFile) {
	*rParamsFile << "\t\t\t<PatchSizeLimit>"
			<< mPatchSizeLimit << "</PatchSizeLimit>\n";
	// Call method on direct parent class
	AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(
			rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ChondrocytesOnlyCellCycleModel)
