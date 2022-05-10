#ifndef PSPC_ITERATOR_TEST_H
#define PSPC_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/iterator/SISIterator.h>
#include <pspc/field/BFieldComparison.h>
#include <util/format/Dbl.h>

#include <fstream>
#include <sstream>


using namespace Util;
using namespace Pscf;
using namespace Pspc;

class SISIteratorTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  setVerbose(0); }

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

   void testConstructors()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      SISIterator<3> sisIter(system);

   }

   template <int D>
   void SetUpSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate1D_lam_rigid.log");

      System<1> system;
      SetUpSystem<1>(system, "in/lam/param.rigid");

      // Read w fields
      system.readWBasis("in/lam/omega.ref");

      // Make a copy of the original field
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate1D_lam_rigid_w.bf");
      system.writeCBasis("out/testIterate1D_lam_rigid_c.bf");

      // Compare solution to original fields
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

      // Check stress value
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT(std::abs(stress) < 1.0E-8);

   }

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate2D_hex_rigid.log");

      System<2> system;
      SetUpSystem<2>(system, "in/hex/param.rigid");

      // Read reference solution
      system.readWBasis("in/hex/omega.ref");

      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Iterate, output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate2D_hex_rigid_w.bf");
      system.writeCBasis("out/testIterate2D_hex_rigid_c.bf");

      // Compare current solution to reference solution
      BFieldComparison comparison(1);
      comparison.compare(wFields_check, system.wFields());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);
      // Maximum error of 2.608E-7 occurs for the first star

      // Check stress
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT (std::abs(stress) < 1.0E-8);
   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testIterate3D_bcc_rigid.log");

      System<3> system;
      SetUpSystem<3>(system, "in/bcc/param.rigid");

      // Read initial guess
      system.readWBasis("in/bcc/omega.ref");

      // Save copy of initial fields
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();

      // Iterate and output solution
      int error = system.iterate();
      if (error) {
         TEST_THROW("Iterator failed to converge.");
      }
      system.writeWBasis("out/testIterate3D_bcc_rigid_w.bf");
      system.writeCBasis("out/testIterate3D_bcc_rigid_c.bf");

      // Compare solution to reference solution
      BFieldComparison comparison(1); // Constructor argument 1 skips star 0
      comparison.compare(wFields_check, system.wFields());
      // setVerbose(1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << comparison.maxDiff() << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 5.0E-7);

      // Test that stress is small
      system.mixture().computeStress();
      double stress = system.mixture().stress(0);
      if (verbose() > 0) {
         std::cout << "stress = " << stress << "\n";
      }
      TEST_ASSERT(std::abs(stress) < 1.0E-7);

   }

};



TEST_BEGIN(SISIteratorTest)
TEST_ADD(SISIteratorTest, testConstructors)
TEST_ADD(SISIteratorTest, testIterate1D_lam_rigid)
TEST_ADD(SISIteratorTest, testIterate2D_hex_rigid)
TEST_ADD(SISIteratorTest, testIterate3D_bcc_rigid)
TEST_END(SISIteratorTest)

#endif
