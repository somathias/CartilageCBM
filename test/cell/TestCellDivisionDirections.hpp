#ifndef TESTCELLDIVISIONDIRECTIONS_HPP_
#define TESTCELLDIVISIONDIRECTIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellDivisionDirection.hpp"
#include "UpwardsCellDivisionDirection.hpp"
#include "DownwardsCellDivisionDirection.hpp"
#include "CellPropertyRegistry.hpp"
#include "TransitCellProliferativeType.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"


class TestCellDivisionDirections: public AbstractCellBasedTestSuite {
public:

	void TestCellDivisionDirectionMethods()  {
		//Test increment/decrement
		MAKE_PTR(UpwardsCellDivisionDirection<3>, p_direction);
		TS_ASSERT_EQUALS(p_direction->GetCellCount(), 0u);

		p_direction->IncrementCellCount();
		TS_ASSERT_EQUALS(p_direction->GetCellCount(), 1u);

		p_direction->DecrementCellCount();
		TS_ASSERT_EQUALS(p_direction->GetCellCount(), 0u);

		TS_ASSERT_THROWS_THIS(p_direction->DecrementCellCount(),
				"Cannot decrement cell count: no cells have this cell property");

		//Test IsType
		TS_ASSERT_EQUALS(p_direction->IsType<TransitCellProliferativeType>(), false);
		TS_ASSERT_EQUALS(p_direction->IsType<UpwardsCellDivisionDirection<3> >(),true);



        MAKE_PTR(UpwardsCellDivisionDirection<3>, p_up_direction);
        MAKE_PTR(DownwardsCellDivisionDirection<3>, p_down_direction);

        //Test IsSame
        TS_ASSERT(p_up_direction->IsSame(p_direction.get()));
        TS_ASSERT(p_direction->IsSame(p_up_direction));

        TS_ASSERT_EQUALS(p_up_direction->IsSame(p_down_direction.get()), false);
        TS_ASSERT_EQUALS(p_down_direction->IsSame(p_up_direction), false);

        // Check that const-ness doesn't matter
//        TS_ASSERT(p_up_direction->IsType<const UpwardsCellDivisionDirection<3>>());
        const UpwardsCellDivisionDirection<3> const_up_direction;

        TS_ASSERT(p_up_direction->IsSame(&const_up_direction));
        TS_ASSERT(const_up_direction.IsSame(p_up_direction));
        TS_ASSERT(const_up_direction.IsSame(p_up_direction.get()));

        //Test colours
        TS_ASSERT_EQUALS(p_up_direction->GetColour(), 17u);
        TS_ASSERT_EQUALS(p_down_direction->GetColour(), 18u);

        //Test directions
        c_vector<double, 3> direction_up = p_up_direction->GetCellDivisionDirection();
        TS_ASSERT_EQUALS(direction_up(0), 0);
        TS_ASSERT_EQUALS(direction_up(1), 0);
        TS_ASSERT_EQUALS(direction_up(2), 1.0);

        c_vector<double, 3> direction_down = p_down_direction->GetCellDivisionDirection();
        TS_ASSERT_EQUALS(direction_down(0), 0);
        TS_ASSERT_EQUALS(direction_down(1), 0);
        TS_ASSERT_EQUALS(direction_down(2), -1.0);


    }

    void TestArchiveUpwardsCellDivisionDirection() 
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "UpwardsCellDivisionDirection.arch";

        // Archive a cell division direction
        {
        	UpwardsCellDivisionDirection<3>* p_direction = new UpwardsCellDivisionDirection<3>();
            p_direction->IncrementCellCount();

            TS_ASSERT_EQUALS(p_direction->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_direction->GetColour(), 17u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_direction;
            output_arch << p_const_state;

            delete p_direction;
        }

        // Restore cell division direction
        {
            AbstractCellProperty* p_direction;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_direction;

            TS_ASSERT_EQUALS(p_direction->GetCellCount(), 1u);

            UpwardsCellDivisionDirection<3>* p_real_state = dynamic_cast<UpwardsCellDivisionDirection<3>*>(p_direction);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 17u);

            // Tidy up
            delete p_direction;
        }
    }

    void TestArchiveDownwardsCellDivisionDirection() 
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DownwardsCellDivisionDirection.arch";

        // Archive a cell division direction
        {
        	DownwardsCellDivisionDirection<3>* p_direction = new DownwardsCellDivisionDirection<3>();
            p_direction->IncrementCellCount();

            TS_ASSERT_EQUALS(p_direction->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_direction->GetColour(), 18u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_direction;
            output_arch << p_const_state;

            delete p_direction;
        }

        // Restore cell division direction
        {
            AbstractCellProperty* p_direction;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_direction;

            TS_ASSERT_EQUALS(p_direction->GetCellCount(), 1u);

            DownwardsCellDivisionDirection<3>* p_real_state = dynamic_cast<DownwardsCellDivisionDirection<3>*>(p_direction);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 18u);

            // Tidy up
            delete p_direction;
        }
    }


};

#endif /*TESTCELLDIVISIONDIRECTIONS_HPP_*/

