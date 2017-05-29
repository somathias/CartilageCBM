/*
 * CellTissueTypeBasedCellCycleModel.hpp
 *
 *  Created on: May 29, 2017
 *      Author: Sonja Mathias
 */

#ifndef CELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_
#define CELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_

#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellTissueType.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"

class CellTissueTypeBasedCellCycleModel: public StochasticDurationGenerationBasedCellCycleModel {
private:

	friend class TestCellTissueTypeBasedCellCycleModel;

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
		archive
				& boost::serialization::base_object<
						StochasticDurationGenerationBasedCellCycleModel>(*this);
	}
public:
	/**
	 * Constructor
	 */
	CellTissueTypeBasedCellCycleModel();

	/**
	 * Overridden builder method to create new copies of
	 * this cell-cycle model.
	 *
	 * @return new cell-cycle model
	 */
	AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set the new cell's tissue type once it has been created after division.
     */
    void InitialiseDaughterCell();

	/**
	 * Outputs cell cycle model parameters to file.
	 *
	 * @param rParamsFile the file stream to which the parameters are output
	 */
	virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellTissueTypeBasedCellCycleModel)

#endif /* CELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_ */
