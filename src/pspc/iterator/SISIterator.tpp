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

   template <int D>
   void SISIterator<D>::setup()
   {

   }

   template <int D>
   int SISIterator<D>::solve()
   {
      // Evaluate scattering functions for this reciprocal space grid
      evaluateScatteringFnc();

      // Solve MDE for initial state and get the W fields 
      system().compute();
      Wfield_ = getWFields();

      // Iterative loop 
      bool done;
      for (int itr = 0; itr < maxItr_; ++itr) {
         
         // Compute error and test it with an isConverged function

         // Do stuff if error is low enough
         if (done) {

            return 0;

         } else {
            // NOTE: This next block of code will likely, in a generalized algorithm,
            // be inside of a loop over the number of monomers (number of W fields).

            // Find the functional derivative with respect to W+ 
            partialPlus_ = findPartialPlus();
            
            // Solve the first semi-implicit equation for W+ (j+1/2)
            WfieldUpdate_[0] = stepWPlus();

            // Update system with these updated fields and re-solve MDEs
            updateSystemFields();
            system().compute();

            // Find the functional derivative with respect to W-.
            partialMinus_ = findPartialMinus();

            // Solve the second semi-implicit equation for W- (j+1/2)
            WfieldUpdate_[1] = stepWMinus();

            // Complete the full step by shifting the spatial average to zero
            shiftAverageZero(WFieldUpdate_)

            // Update system and solve MDEs
            updateSystemFields();
            system().compute();
         }


      }

      // Failure: iteration counter itr reached maxItr without converging
      Log::file() << "Iterator failed to converge before the maximum number of iterations.\n";
      return 1


   }

   template <int D>
   void SISIterator<D>::evaluateScatteringFnc()
   {
      // Evalute diblock scattering functions in Fourier space. This 
      // should only need to be done once as they are only functions of 
      // Fourier variable k, though this may change in a variable unit cell

        
   }

   template <int D>
   void SISIterator<D>::shiftAverageZero(DArray<FieldCPU> & fields)
   {

   }

   template <int D>
   FieldCPU SISIterator<D>::findPartialPlus()
   {

   }

   template <int D>
   FieldCPU SISIterator<D>::findPartialMinus()
   {

   }

   template <int D>
   FieldCPU SISIterator<D>::stepWPlus()
   {

   }

   template <int D>
   FieldCPU SISIterator<D>::stepWMinus()
   {

   }

   template <int D>
   void SISIterator<D>::updateSystemFields()
   {
        
   }

   

}
}
#endif
