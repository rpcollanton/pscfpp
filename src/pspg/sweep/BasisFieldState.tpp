#ifndef PSPG_BASIS_FIELD_STATE_TPP
#define PSPG_BASIS_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BasisFieldState.h"
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
   BasisFieldState<D>::BasisFieldState()
    : FieldState<D, RDField<D> >()
   {}
  
   /*
   * Constructor.
   */
   template <int D>
   BasisFieldState<D>::BasisFieldState(System<D>& system)
    : FieldState<D, RDField<D> >(system)
   {}

   /*
   * Destructor.
   */
   template <int D>
   BasisFieldState<D>::~BasisFieldState()
   {}

   /*
   * Allocate all fields.
   */
   template <int D>
   void BasisFieldState<D>::allocate()
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

      int nMesh = system().mesh().size();
      UTIL_CHECK(nMesh > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (field(i).isAllocated()) {
            UTIL_CHECK(field(i).capacity() == nMesh);
         } else {
            field(i).allocate(nMesh);
         }
      }
   }
 
   /**
   * Read fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::read(const std::string & filename)
   {
      allocate();
      system().fieldIo().readFieldsBasis(filename, fields());
   }

   /**
   * Write fields in symmetry-adapted basis format. 
   */
   template <int D>
   void BasisFieldState<D>::write(const std::string & filename)
   {
      system().fieldIo().writeFieldsBasis(filename, fields(), unitCell());
   }

   /*
   * Get current state of associated System.
   */
   template <int D>
   void BasisFieldState<D>::getSystemState()
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

      int i, j;
      for (i = 0; i < nMonomer; ++i) {
         RDField<D>& stateField = field(i);
         const RDField<D>& systemField = system().wFieldRGrid(i);
         assignReal<<<nBlocks, nThreads>>>
               (stateField.cDField(), systemField.cDField(), nMesh);
      }

   }

   /*
   * Set System state to current state of the BasisFieldState object.
   */
   template <int D>
   void BasisFieldState<D>::setSystemState(bool isFlexible)
   {
      system().setWRGrid(fields());

      if (isFlexible) {
         system().setUnitCell(unitCell());
      }

   }

} // namespace Pspg
} // namespace Pscf
#endif
