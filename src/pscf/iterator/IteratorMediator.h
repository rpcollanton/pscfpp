#ifndef PSCF_ITERATOR_MEDIATOR_H
#define PSCF_ITERATOR_MEDIATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <pscf/AbstractSystem.h>

namespace Pscf {

   // Forward declarations
   class AbstractSystem;
   template <typename T> class Iterator;

   using namespace Util;

   template <typename T>
   class IteratorMediator 
   {
   public:

      /// Constructor
      IteratorMediator(AbstractSystem& sys, Iterator<T>& iter);

      /// Destructor
      ~IteratorMediator(); 

      /// Checks if the system has an initial guess
      virtual bool hasInitialGuess() = 0;
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated 
      virtual int nElements() = 0;

      /// Gets a reference to the current state of the system
      virtual T& getCurrent() = 0;

      /// Runs calculation to evaluate function for fixed point.
      virtual void evaluate() = 0;

      /// Gets residual values from system
      virtual T getResidual() = 0;

      /// Updates the system with a passed in state of the iterator.
      virtual void update(T& newField) = 0;


   private:

      // Const pointer to non-const system
      AbstractSystem * const sys_;

      // Const pointer to non-const iterator
      Iterator * const iter_;

   };

}
#endif