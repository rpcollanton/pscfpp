#ifndef PSPC_ITERATOR_TPP
#define PSPC_ITERATOR_TPP

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

   using namespace Util;

   template<int D>
   Iterator<D>::Iterator(System<D>& system)
    : isFlexible_(false),
      systemPtr_(&system)
   {  setClassName("Iterator"); }

   template<int D>
   Iterator<D>::~Iterator()
   {}

} // namespace Pspc
} // namespace Pscf
#endif
