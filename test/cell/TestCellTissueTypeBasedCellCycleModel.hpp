#ifndef TESTCELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_
#define TTESTCELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "CellTissueTypeBasedCellCycleModel.hpp"
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

class TestCellTissueTypeBasedCellCycleModel: public AbstractCellBasedTestSuite {
public:

	void TestProgressionThroughCellCycle()  {
		TS_ASSERT_THROWS_NOTHING(CellTissueTypeBasedCellCycleModel cell_model);

		unsigned num_cells = (unsigned) 1e5;
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);
		MAKE_PTR(LowerPerichondrialLayer, p_lower_layer);
		MAKE_PTR(HorizontalCellDivisionDirection<3>, p_horizontal);
		MAKE_PTR(UpwardsCellDivisionDirection<3>, p_upwards);

        for (unsigned i=0; i<num_cells; i++)
        {
            CellTissueTypeBasedCellCycleModel* p_cell_cycle_model = new CellTissueTypeBasedCellCycleModel;
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

		double expected_mean_g1_duration = static_cast<CellTissueTypeBasedCellCycleModel*>(cells[0]->GetCellCycleModel())->GetStemCellG1Duration();
        double sample_mean_g1_duration = 0.0;

        for (unsigned i=0; i<num_cells; i++)
        {
            sample_mean_g1_duration += static_cast<CellTissueTypeBasedCellCycleModel*>(cells[i]->GetCellCycleModel())->GetG1Duration()/ (double) num_cells;
        }

        TS_ASSERT_DELTA(sample_mean_g1_duration, expected_mean_g1_duration, 0.1);

		RandomNumberGenerator::Instance()->Reseed(0);
		CellTissueTypeBasedCellCycleModel* p_my_model =
				new CellTissueTypeBasedCellCycleModel;
		// Change G1 Duration for this model
		p_my_model->SetStemCellG1Duration(1.0);

		CellPtr p_my_cell(new Cell(p_state, p_my_model));
		p_my_cell->SetCellProliferativeType(p_stem_type);
		p_my_cell->GetCellData()->SetItem("patch size", 1);
		p_my_cell->AddCellProperty(p_perichondrial_type);
		p_my_cell->AddCellProperty(p_lower_layer);
		p_my_cell->AddCellProperty(p_horizontal);
		p_my_cell->InitialiseCellCycleModel();


		// TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
		// 									PerichondrialCellTissueType>()->GetCellCount(), 1u);
		// TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
		// 									StemCellProliferativeType>()->GetCellCount(), 1u);

		unsigned num_steps = 100;
		double mean_cell_cycle_time =
				p_my_model->GetTransitCellG1Duration() + p_my_model->GetSG2MDuration();

		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(
				mean_cell_cycle_time, num_steps);

		for (unsigned i = 0; i < num_steps; i++) {
			SimulationTime::Instance()->IncrementTimeOneStep();

			//CheckReadyToDivideAndPhaseIsUpdated(p_my_model, 2.35762);
			//std::cout << "Waiting to divide" <<std::endl;

			if (p_my_cell->ReadyToDivide()) {

				std::cout << "Ready to divide" <<std::endl;
				TS_ASSERT_EQUALS(p_perichondrial_type->GetCellCount(), 1u);

				TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>()->GetCellCount(), 0u);


				CellPtr p_daughter_cell = p_my_cell->Divide();
				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								PerichondrialCellTissueType>(), false);
				TS_ASSERT_EQUALS(
						p_my_cell->HasCellProperty<
								PerichondrialCellTissueType>(), true);
				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								ChondrocyteCellTissueType>(), true);
				TS_ASSERT_EQUALS(
						p_my_cell->HasCellProperty<
								ChondrocyteCellTissueType>(), false);

				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								StemCellProliferativeType>(), false);
				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								TransitCellProliferativeType>(), true);
				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								HorizontalCellDivisionDirection<3>>(), false);
				TS_ASSERT_EQUALS(
						p_daughter_cell->HasCellProperty<
								UpwardsCellDivisionDirection<3>>(), true);

				TS_ASSERT_EQUALS(p_perichondrial_type->GetCellCount(), 1u);
//				TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
//													PerichondrialCellTissueType>()->GetCellCount(), 1u);
				TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>()->GetCellCount(), 1u);
			}

		}

	}

	void TestPatchSizeLimit(){
		MAKE_PTR(WildTypeCellMutationState, p_state);

		CellTissueTypeBasedCellCycleModel* p_cell_cycle_model = new CellTissueTypeBasedCellCycleModel;

		TS_ASSERT_EQUALS(p_cell_cycle_model->GetPatchSizeLimit(), 6);

		p_cell_cycle_model->SetPatchSizeLimit(4);

		TS_ASSERT_EQUALS(p_cell_cycle_model->GetPatchSizeLimit(), 4);

		CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

		TS_ASSERT_EQUALS(static_cast<CellTissueTypeBasedCellCycleModel*>(p_cell->GetCellCycleModel())->GetPatchSizeLimit(), 4);

	}


	void TestArchiveCellTissueTypeBasedCellCycleModel()  {

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);
		MAKE_PTR(ChondrocyteCellTissueType, p_chondrocyte_type);

		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath()
				+ "CellTissueTypeBasedCellCycleModel.arch";

		{
			SimulationTime::Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);
			SimulationTime* p_simulation_time = SimulationTime::Instance();
			p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 4);

			CellTissueTypeBasedCellCycleModel* p_model = new CellTissueTypeBasedCellCycleModel;
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

#endif /*TESTCELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_*/
