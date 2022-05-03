#ifndef PSPC_SIS_ITERATOR_TPP
#define PSPC_SIS_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <math.h> // for pow
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
      // Get relevant parameters from system
      const int nField = system().mixture().nMonomer();
      const int nPolymer = system().mixture().nPolymer();
      const int nSolvent = system().mixture().nSolvent(); 
      const IntVec<D> meshDim = system().mesh().dimensions();

      // Make sure SIS is not being used on an inappropriate system!
      if (nField!=2) {
         UTIL_THROW("SIS is currently built only to handle neat 
                        single-component diblock copolymer melts")
      }
      
      // Allocate real-space arrays. Note that wField_[0] does not
      // correspond with the chemical potential field of monomer 0, and is 
      // rather W+ = 1/2*(W_A + W_B). wField_[1] = W- = 1/2*(W_A - W_B)
      wField_.allocate(nField);
      wFieldUpdate_.allocate(nField);
      for (int i = 0; i < nField; i++) {
         wField_[i].allocate(meshDim)
         wFieldUpdate_[i].allocate(meshDim);
      }
      partialPlus_.allocate(meshDim);
      partialMinus_.allocate(meshDim);
      
      // Allocate k-space arrays
      gAA_.allocate(meshDim);
      gAB_.allocate(meshDim);
      gBB_.allocate(meshDim);

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
         // COMPUTE ERROR FUNCTION

         // Do stuff if error is low enough
         if (done) {

            return 0;

         } else {
            // NOTE: This next block of code will likely, in a generalized algorithm,
            // be inside of a loop over the number of monomers (number of W fields).
            // Each step would use the (j+1/2) of all fields previous and (j) of the
            // current field, in a Gauss-Siedel fashion. We would need a generalized
            // "findPartial" function that depends only on the coefficient matrix to 
            // transfrom from one set of fields to another

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

      const IntVec<D> dftDim = gAA_.dftDimensions();
      const FArray<double,6> cellParam = system().unitCell().parameters()
      const double f = system().mixture().polymer(0).block(0).length() / 
                       system().mixture().polymer(0).length();
      
      // Determine actual dft mesh size, accounting for it being cut in half
      int size = 1;
      for (int d = 0; d < D; d++) {
         dftSize*=dftDim[d];
      }
      
      for (int i = 0; i < dftSize i++) {
         temp = i;
         int sizeSlice = dftSize/dftDim[D-1]
         IntVec<D> coord;
         coord[D-1] = size % sizeSlice;
         
         // For each dimension, determine coordinate in k-space array
         for (d = D-2; d >= 0; d++) {
            temp -= coord[d+1]*sizeSlice
            sizeSlice = dftSize/dftDim[d];
            coord[d] = temp % sizeSlice;
         }
         
         // find k norm by summing up over dimensions
         double ksq = 0;
         for (int d = 0; d < D; d++) {
            ksq += pow((double)coord[d]/cellParam[d],2);
         }
         double k = pow(ksq,0.5);
         

         gAA_[i] = 2/pow(k,4) * (f*pow(k,2) + exp(-pow(k,2)*f) - 1);
         gAB_[i] = 1/pow(k,4) * (1 - exp(-pow(k,2)*f))*(1 - exp(-pow(k,2)*(1-f)));
         gBB_[i] = 2/pow(k,4) * ( (1-f)*pow(k,2) + exp(-pow(k,2)*(1-f)) - 1)
      }
      
   }

   template <int D>
   void SISIterator<D>::shiftAverageZero(DArray<RField<D>> & fields)
   {
      // number of fields
      nFields = fields.capacity();

      // for each field
      for (int i = 0; i < nFields; i++) {
         double avg = 0;
         n = fields[i].capacity();
         
         // find average for this field
         for (int j = 0; j < n; j++) {
            avg += fields[i][j]/n;
         }

         // subtract average from field
         for (int j = 0; j < n; j++) {
            fields[i][j] += fields[i][j]/n;
         }
      }

      return;
   }

   template <int D>
   RField<D> SISIterator<D>::findPartialPlus()
   {

   }

   template <int D>
   RField<D> SISIterator<D>::findPartialMinus()
   {

   }

   template <int D>
   RField<D> SISIterator<D>::stepWPlus()
   {
      // do FFT, solve in fourier space, FFT-inverse back
   }

   template <int D>
   RField<D> SISIterator<D>::stepWMinus()
   {
      // solve algebraically (FFT-inverse of g fncs...)

   }

   template <int D>
   DArray<RField<D>> SISIterator<D>::getWFields()
   {
      // get W fields, add/subtract to get W+ and W-, convert
      // to long vector format (can do we do this? should be able to)
   }

   template <int D>
   void SISIterator<D>::updateSystemFields()
   {
      
   }
   

}
}
#endif
