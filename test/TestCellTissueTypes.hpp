#ifndef TESTCELLTISSUETYPES_HPP_
#define TESTCELLTISSUETYPES_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractCellTissueType.hpp"
#include "PerichondrialCellTissueType.hpp"
#include "ChondrocyteCellTissueType.hpp"
#include "CellPropertyRegistry.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCellTissueTypes : public AbstractCellBasedTestSuite
{
public:

    void TestCellTissueTypeMethods() throw(Exception)
    {
        MAKE_PTR(PerichondrialCellTissueType, p_type);
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        p_type->IncrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

        p_type->DecrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        TS_ASSERT_THROWS_THIS(p_type->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");
        TS_ASSERT_EQUALS(p_type->GetColour(), 23u);

        TS_ASSERT_EQUALS(p_type->IsType<PerichondrialCellTissueType>(), true);
        TS_ASSERT_EQUALS(p_type->IsType<ChondrocyteCellTissueType>(), false);

        MAKE_PTR(PerichondrialCellTissueType, p_perichondrial_type);
        MAKE_PTR(ChondrocyteCellTissueType, p_chondrocyte_type);

        TS_ASSERT(p_perichondrial_type->IsSame(p_type.get()));
        TS_ASSERT(p_type->IsSame(p_perichondrial_type));

        TS_ASSERT_EQUALS(p_perichondrial_type->IsSame(p_chondrocyte_type.get()), false);
        TS_ASSERT_EQUALS(p_chondrocyte_type->IsSame(p_perichondrial_type), false);

        // Check that const-ness doesn't matter
        TS_ASSERT(p_perichondrial_type->IsType<const PerichondrialCellTissueType>());
        const PerichondrialCellTissueType const_stem_type;

        TS_ASSERT(p_perichondrial_type->IsSame(&const_stem_type));
        TS_ASSERT(const_stem_type.IsSame(p_perichondrial_type));
        TS_ASSERT(const_stem_type.IsSame(p_perichondrial_type.get()));
    }

    void TestArchivePerichondrialCellTissueType() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "PerichondrialCellTissueType.arch";

        // Archive a cell tissue type
        {
        	PerichondrialCellTissueType* p_type = new PerichondrialCellTissueType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 23u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell tissue type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            PerichondrialCellTissueType* p_real_state = dynamic_cast<PerichondrialCellTissueType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 23u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveChondrocyteCellTissueType() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ChondrocyteCellTissueType.arch";

        // Archive a cell tissue type
        {
        	ChondrocyteCellTissueType* p_type = new ChondrocyteCellTissueType();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 24u);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Write the cell to the archive
            const AbstractCellProperty* const p_const_state = p_type;
            output_arch << p_const_state;

            delete p_type;
        }

        // Restore cell tissue type
        {
            AbstractCellProperty* p_type;

            // Restore the mutation state
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_type;

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

            ChondrocyteCellTissueType* p_real_state = dynamic_cast<ChondrocyteCellTissueType*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 24u);

            // Tidy up
            delete p_type;
        }
    }


};

#endif /*TESTCELLTISSUETYPES_HPP_*/
