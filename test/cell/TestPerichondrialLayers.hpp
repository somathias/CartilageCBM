#ifndef TESTPERICHONDRIALLAYERS_HPP_
#define TESTPERICHONDRIALLAYERS_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "AbstractPerichondrialLayer.hpp"
#include "UpperPerichondrialLayer.hpp"
#include "LowerPerichondrialLayer.hpp"
#include "CellPropertyRegistry.hpp"

#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestPerichondrialLayers : public AbstractCellBasedTestSuite
{
public:

    void TestPerichondrialLayerSetup() 
    {
        MAKE_PTR(LowerPerichondrialLayer, p_type);
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        p_type->IncrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);

        p_type->DecrementCellCount();
        TS_ASSERT_EQUALS(p_type->GetCellCount(), 0u);

        TS_ASSERT_THROWS_THIS(p_type->DecrementCellCount(),
                "Cannot decrement cell count: no cells have this cell property");
        TS_ASSERT_EQUALS(p_type->GetColour(), 45u);

        TS_ASSERT_EQUALS(p_type->IsType<LowerPerichondrialLayer>(), true);
        TS_ASSERT_EQUALS(p_type->IsType<UpperPerichondrialLayer>(), false);

        MAKE_PTR(LowerPerichondrialLayer, p_lower);
        MAKE_PTR(UpperPerichondrialLayer, p_upper);

        TS_ASSERT(p_lower->IsSame(p_type.get()));
        TS_ASSERT(p_type->IsSame(p_lower));

        TS_ASSERT_EQUALS(p_lower->IsSame(p_upper.get()), false);
        TS_ASSERT_EQUALS(p_upper->IsSame(p_lower), false);

        // Check that const-ness doesn't matter
        TS_ASSERT(p_lower->IsType<const LowerPerichondrialLayer>());
        const LowerPerichondrialLayer const_stem_type;

        TS_ASSERT(p_lower->IsSame(&const_stem_type));
        TS_ASSERT(const_stem_type.IsSame(p_lower));
        TS_ASSERT(const_stem_type.IsSame(p_lower.get()));
    }

    void TestArchiveLowerPerichondrialLayer() 
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "LowerPerichondrialLayer.arch";

        // Archive a cell tissue type
        {
        	LowerPerichondrialLayer* p_type = new LowerPerichondrialLayer();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 45u);

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

            LowerPerichondrialLayer* p_real_state = dynamic_cast<LowerPerichondrialLayer*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 45u);

            // Tidy up
            delete p_type;
        }
    }

    void TestArchiveUpperPerichondrialLayer() 
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "UpperPerichondrialLayer.arch";

        // Archive a cell tissue type
        {
        	UpperPerichondrialLayer* p_type = new UpperPerichondrialLayer();
            p_type->IncrementCellCount();

            TS_ASSERT_EQUALS(p_type->GetCellCount(), 1u);
            TS_ASSERT_EQUALS(p_type->GetColour(), 46u);

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

            UpperPerichondrialLayer* p_real_state = dynamic_cast<UpperPerichondrialLayer*>(p_type);
            TS_ASSERT(p_real_state != NULL);
            TS_ASSERT_EQUALS(p_real_state->GetColour(), 46u);

            // Tidy up
            delete p_type;
        }
    }


};

#endif /*TESTPERICHONDRIALLAYERS_HPP_*/
