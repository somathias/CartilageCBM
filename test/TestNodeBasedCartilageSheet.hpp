#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "NodeBasedCartilageSheet.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "FakePetscSetup.hpp"

class TestNodeBasedCartilageSheet: public AbstractCellBasedTestSuite {
public:
	void TestSheet() throw (Exception) {

		// Construct a new cartilage sheet
		NodeBasedCartilageSheet* p_cartilage_sheet =
				new NodeBasedCartilageSheet();

		// set the sheet dimensions
		p_cartilage_sheet->SetCartilageSheetDimensions(4,3,5);
		p_cartilage_sheet->setMaxCoordinatePerturbation(0.1);
		p_cartilage_sheet->UseRandomSeed();
		// generate the nodes
		p_cartilage_sheet->GenerateNodesOnHCPGrid();

		// setup the cell population
		if(!p_cartilage_sheet->isCellPopulationSetup()){
			p_cartilage_sheet->Setup();
		}
		// setup the cell tissue types and cell division directions
		p_cartilage_sheet->InitialiseTissueLayersAndCellDivisionDirections();
		// setup the initial stem cell configuration
		p_cartilage_sheet->InitialiseBulkStemCellConfiguration(2,1);

		// get the cell population
		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population = p_cartilage_sheet->GetCellPopulation();


		//pass it to the simulator
		OffLatticeSimulation<3> simulator(*cell_population);
		simulator.SetOutputDirectory("NodeBasedCartilageSheet");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(10.0);

		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		simulator.AddForce(p_force);

		simulator.Solve();
	}
};
