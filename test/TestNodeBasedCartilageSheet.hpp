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

		NodeBasedCartilageSheet* p_cartilage_sheet =
				new NodeBasedCartilageSheet();

		boost::shared_ptr<NodeBasedCellPopulation<3> > cell_population =
				p_cartilage_sheet->Setup();

		OffLatticeSimulation<3> simulator(*cell_population);
		simulator.SetOutputDirectory("NodeBasedCartilageSheet");
		simulator.SetSamplingTimestepMultiple(12);
		simulator.SetEndTime(10.0);

		MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
		simulator.AddForce(p_force);

		simulator.Solve();
	}
};
