#ifndef PSPG_RGridFieldState_TEST_H
#define PSPG_RGridFieldState_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/sweep/RGridFieldState.h>
#include <pspg/field/RFieldComparison.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class RGridFieldStateTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {}

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      RGridFieldState<3> rfs1(system);
      RGridFieldState<3> rfs2;

   }

   void testRead()
   {
      printMethod(TEST_FUNC);
      
      System<3> system;
      RGridFieldState<3> rfs(system);
   
      // Setup system
      RGridFieldStateTest::SetUpSystem(system);

      // Read in file one way
      system.readWBasis("in/bcc/omega.ref");
      // Read in file another way
      rfs.read("in/bcc/omega.ref");

      // Compare
      RFieldComparison<3> comparison;
      comparison.compare(system.wFields(), rfs.fields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);

   }

   void testWrite()
   {
      // Write tested with a read/write/read/comparison procedure
      printMethod(TEST_FUNC);

      System<3> system;
      RGridFieldState<3> rfs1(system), rfs2(system);

      // Setup system
      RGridFieldStateTest::SetUpSystem(system);

      // read, write, read
      rfs1.read("in/bcc/omega.ref");
      rfs1.write("out/testRGridFieldStateWrite.ref");
      rfs2.read("out/testRGridFieldStateWrite.ref");

      // compare
      RFieldComparison<3> comparison;
      comparison.compare(rfs1.fields(),rfs2.fields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testGetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      RGridFieldState<3> rfs(system);

      // Setup system
      RGridFieldStateTest::SetUpSystem(system);

      // Read in state using system
      system.readWBasis("in/bcc/omega.ref");
      // get it using 
      rfs.getSystemState();

      // compare
      RFieldComparison<3> comparison;
      comparison.compare(rfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);
   }

   void testSetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      RGridFieldState<3> rfs(system);

      // Setup system
      RGridFieldStateTest::SetUpSystem(system);

      // Read in state using 
      rfs.read("in/bcc/omega.ref");
      // set system state
      rfs.setSystemState(true);

      // compare
      RFieldComparison<3> comparison;
      comparison.compare(rfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);
   }

   void testSetSystem()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      RGridFieldState<3> rfs;

      // Setup system
      RGridFieldStateTest::SetUpSystem(system);
      // Invoke setSystem
      rfs.setSystem(system);
   }

   void SetUpSystem(System<3>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/bcc/param.flex", in);
      system.readParam(in);
      in.close();
   }

   template <int D>
   void RDFieldToDField(DArray<DField<cudaReal>> & out, DArray<RDField<D>> const & in)
   {
      // if not allocated, allocate
      int nField = in.capacity();
      int nPoint = in[0].capacity();
      if (!out.isAllocated()) {
         out.allocate(nField);
         for (int i = 0; i < nField; i++) {
            out[i].allocate(nPoint);
         }
      }

      // Copy
      for (int i = 0; i < nField; i++) {
         cudaMemcpy(out[i].cDField(), in[i].cDField(), nPoint*sizeof(cudaReal), cudaMemcpyDeviceToDevice);
      }
   }

};



TEST_BEGIN(RGridFieldStateTest)

TEST_ADD(RGridFieldStateTest, testConstructor)
TEST_ADD(RGridFieldStateTest, testRead)
TEST_ADD(RGridFieldStateTest, testWrite)
TEST_ADD(RGridFieldStateTest, testGetSystemState)
TEST_ADD(RGridFieldStateTest, testSetSystemState)
TEST_ADD(RGridFieldStateTest, testSetSystem)

TEST_END(RGridFieldStateTest)

#endif
