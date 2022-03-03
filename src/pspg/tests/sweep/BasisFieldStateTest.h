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
   
      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in file one way
      system.readWBasis("in/bcc/omega.ref");
      // Read in file another way
      bfs.read("in/bcc/omega.ref");

      // Convert for comparison
      DArray< DField<cudaReal> > d_sysFields, d_bfsFields;
      RDFieldToDField(d_sysFields, system.wFields());
      RDFieldToDField(d_bfsFields, bfs.fields());

      // Compare
      BFieldComparison comparison;
      comparison.compare(d_sysFields, d_bfsFields);
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);

   }

   void testWrite()
   {
      // Write tested with a read/write/read/comparison procedure
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs1(system), bfs2(system);

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // read, write, read
      bfs1.read("in/bcc/omega.ref");
      bfs1.write("out/testBasisFieldStateWrite.ref");
      bfs2.read("out/testBasisFieldStateWrite.ref");

      // Convert for comparison
      DArray< DField<cudaReal> > d_bfs1, d_bfs2;
      RDFieldToDField(d_bfs1, bfs1.fields());
      RDFieldToDField(d_bfs2, bfs2.fields());

      // compare
      BFieldComparison comparison;
      comparison.compare(d_bfs1,d_bfs2);
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testGetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using system
      system.readWBasis("in/bcc/omega.ref");
      // get it using bfs
      bfs.getSystemState();

      // Convert for comparison
      DArray< DField<cudaReal> > d_sysFields, d_bfsFields;
      RDFieldToDField(d_sysFields, system.wFields());
      RDFieldToDField(d_bfsFields, bfs.fields());

      // compare
      BFieldComparison comparison;
      comparison.compare(d_bfsFields,d_sysFields);
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);
   }

   void testSetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using bfs
      bfs.read("in/bcc/omega.ref");
      // set system state
      bfs.setSystemState(true);

      // Convert for comparison
      DArray< DField<cudaReal> > d_sysFields, d_bfsFields;
      RDFieldToDField(d_sysFields, system.wFields());
      RDFieldToDField(d_bfsFields, bfs.fields());

      // compare
      BFieldComparison comparison;
      comparison.compare(d_bfsFields,d_sysFields);
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 1.0e-10);
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



TEST_BEGIN(BasisFieldStateTest)

TEST_ADD(BasisFieldStateTest, testConstructor)
TEST_ADD(BasisFieldStateTest, testRead)
TEST_ADD(BasisFieldStateTest, testWrite)
TEST_ADD(BasisFieldStateTest, testGetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystem)

TEST_END(BasisFieldStateTest)

#endif
