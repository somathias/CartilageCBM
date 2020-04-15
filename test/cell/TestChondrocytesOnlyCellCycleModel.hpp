#ifndef TESTCHONDROCYTESONLYCELLCYCLEMODEL_HPP_
#define TESTCHONDROCYTESONLYCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "ChondrocytesOnlyCellCycleModel.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"

#include "AbstractCellTissueType.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellPropertyRegistry.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <iostream>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestChondrocytesOnlyCellCycleModel: public AbstractCellBasedTestSuite {
public:

	void TestProgressionThroughCellCycle()  {
		TS_ASSERT_THROWS_NOTHING(ChondrocytesOnlyCellCycleModel cell_model);

		unsigned num_cells = (unsigned) 1e5;
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);
        for (unsigned i=0; i<num_cells; i++)
        {
            ChondrocytesOnlyCellCycleModel* p_cell_cycle_model = new ChondrocytesOnlyCellCycleModel;
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
			p_cell->GetCellData()->SetItem("patch size", 1);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

		double expected_mean_g1_duration = static_cast<ChondrocytesOnlyCellCycleModel*>(cells[0]->GetCellCycleModel())->GetStemCellG1Duration();
        double sample_mean_g1_duration = 0.0;

        for (unsigned i=0; i<num_cells; i++)
        {
            sample_mean_g1_duration += static_cast<ChondrocytesOnlyCellCycleModel*>(cells[i]->GetCellCycleModel())->GetG1Duration()/ (double) num_cells;
        }

        TS_ASSERT_DELTA(sample_mean_g1_duration, expected_mean_g1_duration, 0.1);


	}

	void TestPatchSizeLimit(){
		MAKE_PTR(WildTypeCellMutationState, p_state);

		ChondrocytesOnlyCellCycleModel* p_cell_cycle_model = new ChondrocytesOnlyCellCycleModel;

		TS_ASSERT_EQUALS(p_cell_cycle_model->GetPatchSizeLimit(), 6);

		p_cell_cycle_model->SetPatchSizeLimit(4);

		TS_ASSERT_EQUALS(p_cell_cycle_model->GetPatchSizeLimit(), 4);

		CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

		TS_ASSERT_EQUALS(static_cast<ChondrocytesOnlyCellCycleModel*>(p_cell->GetCellCycleModel())->GetPatchSizeLimit(), 4);

	}


	void TestArchiveChondrocytesOnlyCellCycleModel()  {

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);
		MAKE_PTR(ChondrocyteCellTissueType, p_chondrocyte_type);

		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath()
				+ "ChondrocytesOnlyCellCycleModel.arch";

		{
			SimulationTime::Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);
			SimulationTime* p_simulation_time = SimulationTime::Instance();
			p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 4);

			ChondrocytesOnlyCellCycleModel* p_model = new ChondrocytesOnlyCellCycleModel;
			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_transit_type);
			p_cell->AddCellProperty(p_chondrocyte_type);
			p_cell->GetCellData()->SetItem("patch size", 1);
			p_cell->InitialiseCellCycleModel();

			p_simulation_time->IncrementTimeOneStep();
			p_simulation_time->IncrementTimeOneStep();


			p_model->SetBirthTime(-1.0);
			p_model->ReadyToDivide();


			TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);

			CellPtr const p_const_cell = p_cell;

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);
			output_arch << p_const_cell;
		}

		{
			SimulationTime::Destroy();
			SimulationTime* p_simulation_time = SimulationTime::Instance();
			p_simulation_time->SetStartTime(0.0);
			p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

			CellPtr p_cell;

			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_cell;

			AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(p_cell->GetCellCycleModel());

			TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
			TS_ASSERT_DELTA(p_model->GetAge(), 6.0, 1e-12);
			TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);
		}
	}
};

#endif /*TESTCHONDROCYTESONLYCELLCYCLEMODEL_HPP_*/
