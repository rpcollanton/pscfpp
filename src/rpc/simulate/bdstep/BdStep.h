#ifndef RPC_BD_STEP_H
#define RPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class System;
   template <int D> class BdSimulator;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Rpc_Simulate_BdStep_Module
   */
   template <int D>
   class BdStep : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator<D> object
      */
      BdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~BdStep();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Setup before the beginning of each simulation run.
      */
      virtual void setup();

      /**
      * Take a single Brownian dynamics step.
      */
      virtual void step() = 0;
      
      /**
      * Decide whether cc fields need to be saved for move
      * The default implementation is false
      */
      virtual bool needsCc()
      {  return false; }
      
      /**
      * Decide whether dc fields need to be saved for move
      * The default implementation is false
      */
      virtual bool needsDc()
      { return true; }
      
      /**
      * Log output timing results 
      */
      virtual void outputTimers(std::ostream& out);
      
      /**
      * Clear timers 
      */
      virtual void clearTimers();
      
      /**
      * Return converge status of current move.
      */
      bool isConverge() const;

      // Accessor Functions

      /**
      * Output statistics for this move (at the end of simulation)
      */
      virtual void output();

   protected:

      /*
      * Compressor fail to converge. Set isConverge_ to false.
      */
      void failConverge();
      
      /*
      * Compressor converge successfully. Set isConverge_ to true.
      */
      void successConverge();
      
      /**
      * Get parent System object.
      */
      System<D>& system();

      /**
      * Get parent BdSimulator object.
      */
      BdSimulator<D>& simulator();

      /**
      * Get Random number generator of parent System.
      */
      Random& random();

   private:

      /// Compress status.
      bool  isConverge_;
      
      /// Pointer to parent BdSimulator object
      BdSimulator<D>* simulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator
      Random  *randomPtr_;

   };

   // Protected inline methods
   
   /*
   * Return converge status of current move.
   */
   template <int D>
   inline bool BdStep<D>::isConverge() const
   {  return isConverge_; }
   
   /*
   * Compressor fail to converge. Set isConverge_ to false.
   */
   template <int D>
   inline void BdStep<D>::failConverge()
   {  isConverge_ = false; }
   
   /*
   * Compressor converge successfully. Set isConverge_ to true.
   */
   template <int D>
   inline void BdStep<D>::successConverge()
   {  isConverge_ = true; }

   /*
   * Get parent System object.
   */
   template <int D>
   inline System<D>& BdStep<D>::system()
   {  return *systemPtr_; }

   /*
   * Get parent BdSimulator object.
   */
   template <int D>
   inline BdSimulator<D>& BdStep<D>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get Random number generator.
   */
   template <int D>
   inline Random& BdStep<D>::random()
   {  return *randomPtr_; }

   #ifndef RPC_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class BdStep<1>;
   extern template class BdStep<2>;
   extern template class BdStep<3>;
   #endif

}
}
#endif
