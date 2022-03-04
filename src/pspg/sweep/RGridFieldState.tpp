#ifndef PSPG_RGRID_FIELD_STATE_TPP
#define PSPG_RGRID_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RGridFieldState.h"
#include "FieldState.tpp"
#include <util/global.h>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   RGridFieldState<D>::RGridFieldState()
    : FieldState<D, RDField<D> >()
   {}
  
   /*
   * Constructor.
   */
   template <int D>
   RGridFieldState<D>::RGridFieldState(System<D>& system)
    : FieldState<D, RDField<D> >(system)
   {}

   /*
   * Destructor.
   */
   template <int D>
   RGridFieldState<D>::~RGridFieldState()
   {}

   /*
   * Allocate all fields.
   */
   template <int D>
   void RGridFieldState<D>::allocate()
   {
      // Precondition
      UTIL_CHECK(hasSystem());

      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      if (fields().isAllocated()) {
         UTIL_CHECK(fields().capacity() == nMonomer);
      } else {
         fields().allocate(nMonomer);
      }

      IntVec<D> meshdim = system().mesh().dimensions();
      int nMesh = system().mesh().size();
      UTIL_CHECK(nMesh > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (field(i).isAllocated()) {
            UTIL_CHECK(field(i).meshDimensions() == meshdim);
         } else {
            field(i).allocate(meshdim);
         }
      }

      // TEMPORARY. Set up unit cell here. Should be adjusted later.
      unitCell() = system().unitCell();
   }
   /*
   * Get current state of associated System.
   */
   template <int D>
   void RGridFieldState<D>::getSystemState()
   {
      // Get system unit cell
      unitCell() = system().unitCell();
      // Get system wFields
      allocate();
      int nMonomer = system().mixture().nMonomer();
      int nMesh    = system().mesh().size();

      // GPU Resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      for (int i = 0; i < nMonomer; ++i) {
         RDField<D>& stateField = field(i);
         const RDField<D>& systemField = system().wFieldRGrid(i);
         assignReal<<<nBlocks, nThreads>>>
               (stateField.cDField(), systemField.cDField(), nMesh);
      }

   }

   /*
   * Set System state to current state of the RGridFieldState object.
   */
   template <int D>
   void RGridFieldState<D>::setSystemState(bool isFlexible)
   {
      system().setWRGrid(fields());

      if (isFlexible) {
         system().setUnitCell(unitCell());
      }

   }

} // namespace Pspg
} // namespace Pscf
#endif
