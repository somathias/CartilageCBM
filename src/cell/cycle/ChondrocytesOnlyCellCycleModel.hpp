/*
 * ChondrocytesOnlyCellCycleModel.hpp
 *
 *  Created on: May 29, 2017
 *      Author: Sonja Mathias
 */

#ifndef CHONDROCYTESONLYCELLCYCLEMODEL_HPP_
#define CHONDROCYTESONLYCELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellTissueType.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "RandomNumberGenerator.hpp"

class ChondrocytesOnlyCellCycleModel: public AbstractSimpleGenerationalCellCycleModel {
private:

	friend class TestChondrocytesOnlyCellCycleModel;

	/** Needed for serialization. */
	friend class boost::serialization::access;
	/**
	 * Archive the cell-cycle model, never used directly - boost uses this.
	 *
	 * @param archive the archive
	 * @param version the current version of this class
	 */
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version) {
		archive & boost::serialization::base_object<
						AbstractSimpleGenerationalCellCycleModel>(*this);
		RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
		archive & mPatchSizeLimit;

	}

protected:
	/**
	 * Patch size limit. Cells will cease to divide once limit is reached
	 *
	 * Defaults to 6 in the constructor
	 */
	unsigned mPatchSizeLimit;

public:
	/**
	 * Constructor
	 */
	ChondrocytesOnlyCellCycleModel();

	/**
	 * Overridden builder method to create new copies of
	 * this cell-cycle model.
	 *
	 * @return new cell-cycle model
	 */
	AbstractCellCycleModel* CreateCellCycleModel();

	bool ReadyToDivide();

	void SetG1Duration();

	void SetPatchSizeLimit(unsigned);

	unsigned GetPatchSizeLimit();

	/**
	 * Outputs cell cycle model parameters to file.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ChondrocytesOnlyCellCycleModel)

#endif /* CHONDROCYTESONLYCELLCYCLEMODEL_HPP_ */
