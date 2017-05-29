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

	void TestProgressionThroughCellCycle() throw (Exception) {
		TS_ASSERT_THROWS_NOTHING(CellTissueTypeBasedCellCycleModel cell_model);

		MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(TransitCellProliferativeType, p_transit_type);
		MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);

		RandomNumberGenerator::Instance()->Reseed(0);
		CellTissueTypeBasedCellCycleModel* p_my_model =
				new CellTissueTypeBasedCellCycleModel;
		// Change G1 Duration for this model
		p_my_model->SetStemCellG1Duration(1.0);

		CellPtr p_my_cell(new Cell(p_state, p_my_model));
		p_my_cell->SetCellProliferativeType(p_stem_type);
		p_my_cell->AddCellProperty(p_perichondrial_type);
		p_my_cell->InitialiseCellCycleModel();

//		TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
//											PerichondrialCellTissueType>()->GetCellCount(), 1u);
//		TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
//											StemCellProliferativeType>()->GetCellCount(), 1u);

		unsigned num_steps = 100;
		double mean_cell_cycle_time =
				p_my_cell->GetCellCycleModel()->GetAverageStemCellCycleTime();

		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(
				2*mean_cell_cycle_time, num_steps);

		RandomNumberGenerator::Instance()->Reseed(0);

		double g1_duration =
				p_my_cell->GetCellCycleModel()->GetStemCellG1Duration()
						+ 4 * RandomNumberGenerator::Instance()->ranf();
		for (unsigned i = 0; i < num_steps; i++) {
			SimulationTime::Instance()->IncrementTimeOneStep();

			CheckReadyToDivideAndPhaseIsUpdated(p_my_model, g1_duration);
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

				TS_ASSERT_EQUALS(p_perichondrial_type->GetCellCount(), 1u);
//				TS_ASSERT_EQUALS(p_daughter_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
//													PerichondrialCellTissueType>()->GetCellCount(), 1u);
				TS_ASSERT_EQUALS(p_my_cell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<
						ChondrocyteCellTissueType>()->GetCellCount(), 1u);
			}

		}

	}


	void TestArchiveCellTissueTypeBasedCellCycleModel() throw (Exception) {

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

			AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

			TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.0, 1e-12);
			TS_ASSERT_DELTA(p_model->GetAge(), 6.0, 1e-12);
			TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), S_PHASE);
		}
	}
};

#endif /*TESTCELLTISSUETYPEBASEDCELLCYCLEMODEL_HPP_*/
