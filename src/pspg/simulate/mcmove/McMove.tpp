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
   }

   /*
   * Trivial default implementation - always returns false.
   */
   template <int D>
   bool McMove<D>::move()
   { 
      incrementNAttempt();

      // Get current Hamiltonian
      double oldHamiltonian = mcSimulator().mcHamiltonian();
     
      // Save current state 
      mcSimulator().saveMcState();
      
      // Clear both eigen-components of the fields and mcHamiltonian 
      mcSimulator().clearData();

      // Attempt modification
      attemptMove();

      // Call compressor
      system().compressor().compress();

      // Evaluate new Hamiltonian
      mcSimulator().computeWC();
      mcSimulator().computeMcHamiltonian();
      double newHamiltonian = mcSimulator().mcHamiltonian();
     // Log::file() << "newHamiltonian" << newHamiltonian << "\n";
      //Log::file() << "oldHamiltonian" << oldHamiltonian << "\n";
      // Accept or reject move
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
      // Log::file() << "newidealHamiltonian" << mcSimulator().mcIdealHamiltonian() << "\n";
      // Log::file() << "accept" << accept << "\n";
      return accept;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void McMove<D>::output()
   {}

}
}
#endif
