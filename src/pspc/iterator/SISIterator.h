#ifndef PSPC_SIS_ITERATOR_H
#define PSPC_SIS_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pspc/field/RFieldDft.h>
#include <pspc/field/RField.h>

namespace Pscf
{
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspc implementation of the Semi-Implicit Seidel iterator.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class SISIterator : public Iterator<D>
   {

   public:
      /**
      * Constructor.
      *
      * \param system System object associated with this iterator.
      */
      SISIterator(System<D> &system);

      /**
      * Destructor.
      */
      ~SISIterator();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream &in);

      /**
      * Setup and allocate required memory.
      */
      void setup();

      /**
      * Iterate to a solution
      */
      int solve();

   protected:
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;

   private:
      /// Error tolerance
      double epsilon_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      /// Maximum number of iterations to attempt
      int maxItr_;

      /// Time step size
      double dt_;

      /// W fields last used to solve the MDEs. Has two components,
      /// W+ = 1/2*(W_A + W_B) (pressure-like field) and
      /// W- = 1/2*(W_A - W_B) (exchange field)

      DArray<RField<D>> Wfields_;

      /// Fourier transformed scattering function, gAA
      RFieldDft<D> gAA_;

      /// Fourier transformed scattering function, gAB
      RFieldDft<D> gAB_;

      /// Fourier transformed scattering function, gBB
      RFieldDft<D> gBB_;

      /// Real-space partial functional derivative ofs the effective
      /// Hamiltonian with respect to the transformed W fields
      DArray<RField<D>> partialDeriv_;

      /**
      * Compute the diblock scattering functions.
      */
      void evaluateScatteringFnc();

      /**
      * Shift the average of each field in an array of fields to 0.
      *
      * \param fields An array of fields to be shifted.
      */
      void shiftAverageZero(DArray<RField<D>> &fields);

      /**
      * Find the partial functional derivative of the effective
      * Hamiltonian with respect to the W+ field.
      */
      RField<D> findPartialPlus();

      /**
      * Find the partial functional derivative of the effective
      * Hamiltonian with respect to the W- field.
      */
      RField<D> findPartialMinus(const RField<D> &WMinus);

      /**
      * Solve the semi-implicit equation for the next half-step of W+.
      */
      RField<D> stepWPlus(const RField<D> &WPlus, const RField<D> &partialPlus);

      /**
      * Solve the semi-implicit equation for the next half-step of W+.
      */
      RField<D> stepWMinus(const RField<D> &WMinus, const RField<D> &partialMinus);

      /**
      * Get the current W fields on the system, and convert W_A and W_B
      * to W+ and W- (in that order). Output them.
      */
      void getWFields(DArray<RField<D>> &Wfields);

      /**
      * Checks if the system is converged by comparing the self-consistency
      * equations via taking the partial functional derivatives with respect
      * to W+/-.
      */
      bool isConverged();

      /**
      * Update the W fields on the associated system object
      * with the W fields stored in WfieldsUpdate
      */
      void updateSystemFields(const DArray<RField<D>> &WfieldsUpdate);
   };

} // namespace Pspc
} // namespace Pscf
#endif
