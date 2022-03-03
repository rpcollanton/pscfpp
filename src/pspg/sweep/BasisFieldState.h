#ifndef PSPG_BASIS_FIELD_STATE_H
#define PSPG_BASIS_FIELD_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/math/GpuResources.h>
#include <pspg/field/RDField.h>
#include "FieldState.h"
#include <string>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /**
   * FieldState for fields in symmetry-adapted basis format.
   */
   template <int D>
   class BasisFieldState : public FieldState<D, RDField<D> >
   {
   public:

      /**
      * Default constructor.
      */
      BasisFieldState();

      /**
      * Constructor, create association with a parent system.
      *
      * \param system associated parent system
      */
      BasisFieldState(System<D>& system);

      /**
      * Destructor.
      */
      ~BasisFieldState();

      /**
      * Allocate all fields.
      *
      * Precondition: hasSystem() == true
      */
      void allocate();

      /**
      * Read state from file.
      *
      * \param filename name of input w-field file in symmetry-adapted format.
      */
      void read(const std::string & filename);
   
      /**
      * Write state to file.
      *
      * \param filename name of output file, in symmetry-adapated format.
      */
      void write(const std::string & filename);

      /**
      * Copy the current state of the associated system.
      *
      * Copy the fields and the unit cell.
      */
      void getSystemState();

      /**
      * Set the state of the associated system to this state.
      *
      * \param isFlexible update system unit cell iff isFlexible == true.
      */
      void setSystemState(bool isFlexible);

      // Inherited member functions
      using FieldState<D, RDField<D> >::fields;
      using FieldState<D, RDField<D> >::field;
      using FieldState<D, RDField<D> >::unitCell;
      using FieldState<D, RDField<D> >::system;
      using FieldState<D, RDField<D> >::hasSystem;
      using FieldState<D, RDField<D> >::setSystem;

   };

   #ifndef PSPG_BASIS_FIELD_STATE_TPP
   // Suppress implicit instantiation
   extern template class BasisFieldState<1>;
   extern template class BasisFieldState<2>;
   extern template class BasisFieldState<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
#endif
