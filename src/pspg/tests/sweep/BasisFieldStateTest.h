#ifndef PSPG_BASISFIELDSTATE_TEST_H
#define PSPG_BASISFIELDSTATE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/sweep/BasisFieldState.h>
#include <pspg/field/BFieldComparison.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class BasisFieldStateTest : public UnitTest
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
      BasisFieldState<3> bfs1(system);
      BasisFieldState<3> bfs2;

   }

   void testRead()
   {
      printMethod(TEST_FUNC);
      
      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;
   
      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in file one way
      system.readWBasis("in/bcc/omega.ref");
      // Read in file another way
      bfs.read("in/bcc/omega.ref");
      // Compare
      comparison.compare(bfs.fields(), system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);

   }

   void testWrite()
   {
      // Write tested with a read/write/read/comparison procedure
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs1(system), bfs2(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // read, write, read
      bfs1.read("in/bcc/omega.ref");
      bfs1.write("out/testBasisFieldStateWrite.ref");
      bfs2.read("out/testBasisFieldStateWrite.ref");

      // compare
      comparison.compare(bfs1.fields(),bfs2.fields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testGetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using system
      system.readWBasis("in/bcc/omega.ref");
      // get it using bfs
      bfs.getSystemState();
      // compare
      comparison.compare(bfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testSetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using bfs
      bfs.read("in/bcc/omega.ref");
      // set system state
      bfs.setSystemState(true);
      // compare
      comparison.compare(bfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testSetSystem()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);
      // Invoke setSystem
      bfs.setSystem(system);
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

};



TEST_BEGIN(BasisFieldStateTest)

TEST_ADD(BasisFieldStateTest, testConstructor)
TEST_ADD(BasisFieldStateTest, testRead)
TEST_ADD(BasisFieldStateTest, testWrite)
TEST_ADD(BasisFieldStateTest, testGetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystem)

TEST_END(BasisFieldStateTest)

#endif
