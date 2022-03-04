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
      comparison.compare(rfs.fields(),system.wFieldsRGrid());
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

      // Read basis format file into rgrid field
      rfs.allocate();
      readBasisToRGrid("in/bcc/omega.ref", rfs.fields(), system);
      // set system state
      rfs.setSystemState(true);

      // compare
      RFieldComparison<3> comparison;
      comparison.compare(rfs.fields(), system.wFieldsRGrid());
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
   void readBasisToRGrid(std::string filename, DArray<RDField<D>> & rfields, System<D> & system) 
   {
      // Holding area for basis 
      DArray<RDField<D>> tempBasis;
      int nBasis = system.basis().nStar(); 
      
      // Allocate
      tempBasis.allocate(system.mixture().nMonomer());
      for (int i = 0; i < system.mixture().nMonomer(); i++) {
         tempBasis[i].allocate(nBasis);
      }

      // Read in field in basis format 
      system.fieldIo().readFieldsBasis(filename, tempBasis);

      // Convert
      system.fieldIo().convertBasisToRGrid(tempBasis, rfields);

   }

};



TEST_BEGIN(RGridFieldStateTest)

TEST_ADD(RGridFieldStateTest, testConstructor)
TEST_ADD(RGridFieldStateTest, testGetSystemState)
TEST_ADD(RGridFieldStateTest, testSetSystemState)
TEST_ADD(RGridFieldStateTest, testSetSystem)

TEST_END(RGridFieldStateTest)

#endif
