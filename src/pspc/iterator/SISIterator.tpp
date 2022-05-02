#ifndef PSPC_SIS_ITERATOR_TPP
#define PSPC_SIS_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "SISIterator.h"
#include <pspc/System.h>
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   template <int D>
   void SISIterator<D>::readParameters(std::istream& in)
   {
       // convergence criterion type
       // convergence tolerance
       // max number of iterations
   }

   

}
}
#endif
