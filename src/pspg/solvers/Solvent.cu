/*
* PSCF - Polymer Self-Consistent Field Theory 
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.tpp"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   // Explicit instantiation of relevant class instances
   template class Solvent<1>;
   template class Solvent<2>;
   template class Solvent<3>;

}
}
