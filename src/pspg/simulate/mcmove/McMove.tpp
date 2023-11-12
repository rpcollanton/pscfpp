#ifndef PSPG_MC_MOVE_TPP
#define PSPG_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"
#include <pspg/System.h>
#include <pspg/simulate/McSimulator.h>
#include <util/archives/Serializable_includes.h>
#include <pspg/compressor/Compressor.h>

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMove<D>::McMove(McSimulator<D>& mcSimulator) 
    : mcSimulatorPtr_(&mcSimulator),
      systemPtr_(&(mcSimulator.system())),
      randomPtr_(&(mcSimulator.random())),
      nAttempt_(0),
      nAccept_(0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   McMove<D>::~McMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void McMove<D>::readParameters(std::istream &in)
   {}
   
   /*
   * Read the probability from file.
   */
   template <int D>
   void McMove<D>::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }
  
   /*
   * Setup at beginning of loop.
   *
   * Trivial default implementation - initializes counters.
   */
   template <int D>
   void McMove<D>::setup()
   { 
      nAttempt_ = 0;
      nAccept_  = 0;
      clearTimers();
   }

   /*
   * Trivial default implementation - always returns false.
   */
   template <int D>
   bool McMove<D>::move()
   { 
      totalTimer_.start();
      incrementNAttempt();

      // Get current Hamiltonian
      double oldHamiltonian = mcSimulator().hamiltonian();
     
      // Save current state 
      mcSimulator().saveMcState();
      
      // Clear both eigen-components of the fields and hamiltonian 
      mcSimulator().clearData();

      // Attempt modification
      attemptMoveTimer_.start();
      attemptMove();
      attemptMoveTimer_.stop();
      
      // Call compressor
      compressorTimer_.start();
      system().compressor().compress();
      compressorTimer_.stop();
      
      // Compute eigenvector components of the current w fields
      computeWCTimer_.start();
      mcSimulator().computeWC();
      computeWCTimer_.stop();
      
      // Evaluate new Hamiltonian
      computeHamiltonianTimer_.start();
      mcSimulator().computeHamiltonian();
      double newHamiltonian = mcSimulator().hamiltonian();
      computeHamiltonianTimer_.stop();
      // Log::file() << "newHamiltonian" << newHamiltonian << "\n";
      //Log::file() << "oldHamiltonian" << oldHamiltonian << "\n";
      // Accept or reject move
     
      decisionTimer_.start();
      bool accept = false;
      double weight = exp(-(newHamiltonian - oldHamiltonian));
      //Log::file() << "weight" << weight << "\n";
      accept = random().metropolis(weight);
      if (accept) {
          incrementNAccept();
          mcSimulator().clearMcState();
      } else {
          mcSimulator().restoreMcState();
      }
     /// Log::file() << "newFieldHamiltonian" << mcSimulator().mcFieldHamiltonian() << "\n";
      // Log::file() << "newidealHamiltonian" << mcSimulator().idealHamiltonian() << "\n";
      // Log::file() << "accept" << accept << "\n";
      decisionTimer_.stop();
      totalTimer_.stop();
      return accept;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void McMove<D>::output()
   {}
   
   template<int D>
   void McMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "McMove times contributions:\n";
      // Output timing results, if requested.
      double total = totalTimer_.time();
      out << "                          "
          << "Total" << std::setw(17) << "Per Move" << std::setw(14) << "Fraction" << "\n";
      out << "Attempt Move:             "
          << Dbl(attemptMoveTimer_.time(), 9, 3)  << " s,  "
          << Dbl(attemptMoveTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(attemptMoveTimer_.time()/total, 9, 3) << "\n";
      out << "Compressor:               "
          << Dbl(compressorTimer_.time(), 9, 3)  << " s,  "
          << Dbl(compressorTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(compressorTimer_.time()/total, 9, 3) << "\n";
      out << "Compute eigen-components: "
          << Dbl(computeWCTimer_.time(), 9, 3)  << " s,  "
          << Dbl(computeWCTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(computeWCTimer_.time()/total, 9, 3) << "\n";
      out << "Compute Hamiltonian:      "
          << Dbl(computeHamiltonianTimer_.time(), 9, 3)  << " s,  "
          << Dbl(computeHamiltonianTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(computeHamiltonianTimer_.time()/total, 9, 3) << "\n";
      out << "Accept or Reject:         "
          << Dbl(decisionTimer_.time(), 9, 3)  << " s,  "
          << Dbl(computeHamiltonianTimer_.time()/nAttempt_, 9, 3)  << " s,  "
          << Dbl(decisionTimer_.time()/total, 9, 3) << "\n";
      out << "total time:               "
          << Dbl(total, 9, 3) << " s,  "
          << Dbl(total/nAttempt_, 9, 3) << " s  \n";
      out << "\n";
   }
   
   template<int D>
   void McMove<D>::clearTimers()
   {
      attemptMoveTimer_.clear();
      compressorTimer_.clear();
      computeWCTimer_.clear();
      computeHamiltonianTimer_.clear();
      decisionTimer_.clear();
      totalTimer_.clear();
   }
      

}
}
#endif
