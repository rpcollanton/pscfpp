#ifndef PSPC_SIMULATOR_TEST_H
#define PSPC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/simulate/Simulator.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Pspc;

class SimulatorTest : public LogFileUnitTest
{

   System<3> system;

public:


   SimulatorTest()
    : system()
   {}

   void setUp()
   {  setVerbose(0); }

   void initSystem(std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testSystem.log");
      ParamComponent::setEcho(true);

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

   }

   void testAnalyzeChi()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param1_simulator");
      Simulator<3> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      double chi = system.interaction().chi(0,1);
      TEST_ASSERT( fabs(system.interaction().chi(0,0)) < 1.0E-8);
      TEST_ASSERT( fabs(system.interaction().chi(1,1)) < 1.0E-8);

      DArray<double> vals = simulator.chiEvals();
      TEST_ASSERT( fabs((vals[0] - simulator.chiEval(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((vals[1] - simulator.chiEval(1))/chi) < 1.0E-8);
      TEST_ASSERT(fabs((vals[0] + chi)/chi) < 1.0E-8);
      TEST_ASSERT(fabs(vals[1]/chi) < 1.0E-8);

      DMatrix<double> vecs = simulator.chiEvecs();
      TEST_ASSERT(fabs(vecs(0,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(0,1) + 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,1) - 1.0) < 1.0E-8);

      DArray<double>  sc = simulator.sc();
      TEST_ASSERT( fabs((sc[0] - simulator.sc(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((sc[1] - simulator.sc(1))/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(0)/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(1)/chi - 0.5) < 1.0E-8);

      #if 0
      std::cout << std::endl;
      std::cout << "vals  = " << vals[0] << "  " << vals[1] << std::endl;
      std::cout << "vec0  = " << vecs(0,0) << "  " << vecs(0,1) << std::endl;
      std::cout << "vec1  = " << vecs(1,0) << "  " << vecs(1,1) << std::endl;
      #endif

   }

   void testSaddlePointField()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param1_simulator");
      Simulator<3> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      system.readWRGrid("in/w_gyr.rf");
      DArray< RField<3> > const & w = system.w().rgrid();

      system.compute();
      DArray< RField<3> > const & c = system.c().rgrid();

      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();

      simulator.computeWc();
      DArray< RField<3> > const & wc = simulator.wc();

      simulator.computeCc();
      DArray< RField<3> > const & cc = simulator.cc();

      simulator.computeDc();
      DArray< RField<3> > const & dc = simulator.dc();

      // Check allocation and capacities
      TEST_ASSERT(c.capacity() == nMonomer);
      TEST_ASSERT(w.capacity() == nMonomer);
      TEST_ASSERT(wc.capacity() == nMonomer);
      TEST_ASSERT(cc.capacity() == nMonomer);
      int i;
      for (i=0; i < nMonomer; ++i) {
         TEST_ASSERT(c[i].capacity() == meshSize);
         TEST_ASSERT(w[i].capacity() == meshSize);
         TEST_ASSERT(wc[i].capacity() == meshSize);
         TEST_ASSERT(cc[i].capacity() == meshSize);
      }
      TEST_ASSERT(dc.capacity() == nMonomer - 1);
      for (i=0; i < nMonomer - 1; ++i) {
         TEST_ASSERT(dc[i].capacity() == meshSize);
      }

      // Test wc field
      for (i = 0; i < meshSize; ++i) {
         TEST_ASSERT(fabs( w[0][i] - wc[0][i] - wc[1][i] ) < 1.0E-6);
         TEST_ASSERT(fabs( w[0][i] - w[1][i] - 2.0*wc[0][i] ) < 1.0E-6);
      }

      // Test cc field
      for (i = 0; i < meshSize; ++i) {
         //std::cout << c[0][i] - c[1][i] << "   " << cc[0][i] << std::endl;
         TEST_ASSERT(fabs( c[0][i] - c[1][i] - cc[0][i] ) < 1.0E-6);
      }

      // Test dc field
      for (i = 0; i < meshSize; ++i) {
         //std::cout << dc[0][i] << std::endl;
         TEST_ASSERT(fabs( dc[0][i] ) < 1.0E-5);
      }

      #if 0
      std::cout << std::endl;
      std::cout << system.domain().mesh().dimensions() << std::endl;
      std::cout << system.domain().mesh().size() << std::endl;
      #endif

      system.computeFreeEnergy();
      double fHelmholtz = system.fHelmholtz();
      double volume = system.domain().unitCell().volume();
      double vMonomer = system.mixture().vMonomer();
      double ratio = volume/vMonomer;
      simulator.computeHamiltonian();
      double hamiltonian = simulator.hamiltonian();

      TEST_ASSERT( fabs((hamiltonian - ratio*fHelmholtz)/hamiltonian) < 1.0E-6);

      #if 0
      std::cout << "fHelmholtz = " << fHelmholtz << "  " 
                << ratio*fHelmholtz  << std::endl;

      std::cout << "Hamiltonian = " << hamiltonian/ratio << "  " 
                << hamiltonian << std::endl;
      #endif
     
   }

};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testAnalyzeChi)
TEST_ADD(SimulatorTest, testSaddlePointField)
TEST_END(SimulatorTest)

#endif
