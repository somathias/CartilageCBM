#ifndef TESTCELLDIVISIONDIRECTIONSANDTISSUETYPESWRITERS_HPP_
#define TESTCELLDIVISIONDIRECTIONSANDTISSUETYPESWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SimulationTime.hpp"

// Cell writers
#include "CellTissueTypesWriter.hpp"
#include "CellDivisionDirectionsWriter.hpp"

// tissue types
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"

#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"

#include <iostream>

#include "PetscSetupAndFinalize.hpp"
//#include "FakePetscSetup.hpp"

// Note that high level tests of all cell writers can be found in the
// TestMeshBasedCellPopulation::TestMeshBasedCellPopulationWriteResultsToFile

class TestCellDivisionDirectionsAndTissueTypesWriters: public AbstractCellBasedTestSuite {
public:

	void TestCellTissueTypesWriter() throw (Exception) {
		EXIT_IF_PARALLEL;

		// Set up SimulationTime (this is usually done by a simulation object)
		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

		// Create a simple node-based cell population
		std::vector<Node<2>* > nodes;
		nodes.push_back(new Node<2>(0, false, 1.4));
		nodes.push_back(new Node<2>(1, false, 2.3));
		nodes.push_back(new Node<2>(2, false, -6.1));

		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_perichondrial(CellPropertyRegistry::Instance()->Get<PerichondrialCellTissueType>());
		boost::shared_ptr<AbstractCellProperty> p_chondrocyte(CellPropertyRegistry::Instance()->Get<ChondrocyteCellTissueType>());

		std::vector<CellPtr> cells;
		for (unsigned i=0; i<3; i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
			CellPtr p_cell(new Cell(p_state, p_cell_model));
			p_cell->SetBirthTime(-0.7 - i*0.5);
			cells.push_back(p_cell);
		}

		cells[0]->AddCellProperty(p_chondrocyte);
		cells[1]->AddCellProperty(p_perichondrial);

		NodeBasedCellPopulation<2> cell_population(mesh, cells);

		// Create output directory
		std::string output_directory = "TestCellTissueTypesWriter";
		OutputFileHandler output_file_handler(output_directory, false);
		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

		// Create cell writer and output data for each cell to file
		CellTissueTypesWriter<2,2> cell_writer;
		cell_writer.OpenOutputFile(output_file_handler);
		cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator other_cell_iter = cell_population.Begin();
				other_cell_iter != cell_population.End();
				++other_cell_iter)
		{
			cell_writer.VisitCell(*other_cell_iter, &cell_population);
		}
		cell_writer.CloseFile();

		// Test that the data are output correctly
		FileComparison(results_dir + "results.vizcelltissuetypes", "projects/chaste_project/test/data/TestCellDivisionDirectionsAndTissueTypesWriters/results.vizcelltissuetypes").CompareFiles();

		// Test the correct data are returned for VTK output for the first cell
		double vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(0), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 24.0, 1e-6);

		vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(1), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 23.0, 1e-6);

		vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(2), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

		// Test GetVtkCellDataName() method
		TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell tissue types");

		// Avoid memory leak
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void TestCellTissueTypesWriterArchiving() throw (Exception)
	{
		// The purpose of this test is to check that archiving can be done for this class
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellTissueTypesWriter.arch";

		{
			AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellTissueTypesWriter<2,2>();

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			output_arch << p_cell_writer;

			delete p_cell_writer;
		}
		PetscTools::Barrier(); //Processes read after last process has (over-)written archive
		{
			AbstractCellBasedWriter<2,2>* p_cell_writer_2;

			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_cell_writer_2;

			delete p_cell_writer_2;
		}
	}

	void TestCellDivisionDirectionsWriter() throw (Exception) {
		EXIT_IF_PARALLEL;

		// Set up SimulationTime (this is usually done by a simulation object)
		SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

		// Create a simple node-based cell population
		std::vector<Node<2>* > nodes;
		nodes.push_back(new Node<2>(0, false, 1.4));
		nodes.push_back(new Node<2>(1, false, 2.3));
		nodes.push_back(new Node<2>(2, false, -6.1));

		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, 1.5);

		boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_upwards(CellPropertyRegistry::Instance()->Get<UpwardsCellDivisionDirection<2> >());
		boost::shared_ptr<AbstractCellProperty> p_downwards(CellPropertyRegistry::Instance()->Get<DownwardsCellDivisionDirection<2> >());

		std::vector<CellPtr> cells;
		for (unsigned i=0; i<3; i++)
		{
			FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
			CellPtr p_cell(new Cell(p_state, p_cell_model));
			p_cell->SetBirthTime(-0.7 - i*0.5);
			cells.push_back(p_cell);
		}

		cells[0]->AddCellProperty(p_upwards);
		cells[1]->AddCellProperty(p_downwards);

		NodeBasedCellPopulation<2> cell_population(mesh, cells);

		// Create output directory
		std::string output_directory = "TestCellDivisionDirectionsWriter";
		OutputFileHandler output_file_handler(output_directory, false);
		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

		// Create cell writer and output data for each cell to file
		CellDivisionDirectionsWriter<2,2> cell_writer;
		cell_writer.OpenOutputFile(output_file_handler);
		cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator other_cell_iter = cell_population.Begin();
				other_cell_iter != cell_population.End();
				++other_cell_iter)
		{
			cell_writer.VisitCell(*other_cell_iter, &cell_population);
		}
		cell_writer.CloseFile();

		// Test that the data are output correctly
		FileComparison(results_dir + "results.vizcelldivisiondirections", "projects/chaste_project/test/data/TestCellDivisionDirectionsAndTissueTypesWriters/results.vizcelldivisiondirections").CompareFiles();

		// Test the correct data are returned for VTK output for the first cell
		double vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(0), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 17.0, 1e-6);

		vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(1), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 18.0, 1e-6);

		vtk_data = cell_writer.GetCellDataForVtkOutput(cell_population.GetCellUsingLocationIndex(2), &cell_population);
		TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

		// Test GetVtkCellDataName() method
		TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell division directions");

		// Avoid memory leak
		for (unsigned i=0; i<nodes.size(); i++)
		{
			delete nodes[i];
		}
	}

	void TestCellDivisionDirectionsWriterArchiving() throw (Exception)
	{
		// The purpose of this test is to check that archiving can be done for this class
		OutputFileHandler handler("archive", false);
		std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellDivisionDirectionsWriter.arch";

		{
			AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellTissueTypesWriter<2,2>();

			std::ofstream ofs(archive_filename.c_str());
			boost::archive::text_oarchive output_arch(ofs);

			output_arch << p_cell_writer;

			delete p_cell_writer;
		}
		PetscTools::Barrier(); //Processes read after last process has (over-)written archive
		{
			AbstractCellBasedWriter<2,2>* p_cell_writer_2;

			std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
			boost::archive::text_iarchive input_arch(ifs);

			input_arch >> p_cell_writer_2;

			delete p_cell_writer_2;
		}
	}

};

#endif /*TESTCELLDIVISIONDIRECTIONSANDTISSUETYPESWRITERS_HPP_*/
