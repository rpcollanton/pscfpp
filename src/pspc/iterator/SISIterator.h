#ifndef PSPC_SIS_ITERATOR_H
#define PSPC_SIS_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspc implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class SISIterator : Iterator<D>
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this iterator.
      */
      SISIterator(System<D>& system);

      /**
      * Destructor.
      */
      ~SISIterator();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Setup and allocate required memory.
      */
      void setup();

      /**
      * Iterate to a solution
      */
      int solve();

   protected:
   
      using ParamComposite::readOptional;
      using Iterator<D>::system;

   private:

      /// Error tolerance
      double epsilon_;

      /// Maximum number of iterations to attempt
      int maxItr_;

      /// W fields last used to solve the MDEs. Has two components,
      /// W+ = 1/2*(W_A + W_B) and W- = 1/2*(W_A - W_B)
      DArray<FieldCPU> Wfield_;

      /// W fields updated with SIS and average shift
      DArray<FieldCPU> WfieldUpdate_;

      /// Fourier transformed scattering function, gAA
      FieldCPU gAA_;

      /// Fourier transformed scattering function, gAB
      FieldCPU gAB_;

      /// Fourier transformed scattering function, gBB
      FieldCPU gBB_;

      /**
      * Compute the diblock scattering functions.
      */
      void evaluateScatteringFnc();

      /**
      * Shift the average of each field in an array of fields
      * to 0.
      * 
      * \param fields An array of fields to be shifted.
      */
      DArray<FieldCPU> shiftAverageZero(DArray<FieldCPU> fields);

      /**
      * Find the partial functional derivative of the effective
      * Hamiltonian with respect to the W+ field.
      */
      FieldCPU findPartialPlus();

      /**
      * Find the partial functional derivative of the effective
      * Hamiltonian with respect to the W- field.
      */
      FieldCPU findPartialMinus();

   };

} // namespace Pspc
} // namespace Pscf
#endif
