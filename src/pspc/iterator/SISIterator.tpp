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
      // max number of iterations
      read(in, "maxItr", maxItr_);
      // convergence tolerance
      read(in, "epsilon", epsilon_);
      // pseudo-time step size
      read(in, "timeStep", dt_);
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
      if (nField!=2 || nPolymer!=1 || nSolvent!=0) {
         UTIL_THROW("SIS is currently built only to handle neat single-component diblock copolymer melts.");
      }
      
      // Allocate real-space arrays. Note that Wfields_[0] does not
      // correspond with the chemical potential field of monomer 0, and is 
      // rather W+ = 1/2*(W_A + W_B). Wfields_[1] = W- = 1/2*(W_A - W_B)
      Wfields_.allocate(nField);
      WfieldsUpdate_.allocate(nField);
      partialDeriv_.allocate(nField);
      for (int i = 0; i < nField; i++) {
         Wfields_[i].allocate(meshDim);
         WfieldsUpdate_[i].allocate(meshDim);
         partialDeriv_[i].allocate(meshDim);
      }
      
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
      Wfields_ = getWFields();

      // Iterative loop 
      bool done;
      for (int itr = 0; itr < maxItr_; ++itr) {
         
         // Compute error and test it with an isConverged function
         done = isConverged();

         // Do stuff if error is low enough
         if (done) {

            return 0;

         } else {
            // NOTE: This next block of code will likely, in a generalized algorithm,
            // be inside of a loop over the number of monomers (number of W fields).
            // Each step would use the (j+1/2) of all fields previous and (j) of the
            // current field, in a Gauss-Seidel fashion. We would need a generalized
            // "findPartial" function that depends only on the coefficient matrix to 
            // transfrom from one set of fields to another

            // Find the functional derivative with respect to W+ 
            partialDeriv_[0] = findPartialPlus();
            
            // Solve the first semi-implicit equation for W+ (j+1/2)
            WfieldsUpdate_[0] = stepWPlus(Wfields_[0], partialDeriv_[0]);

            // Update system with these updated fields and re-solve MDEs
            updateSystemFields();
            system().compute();

            // Find the functional derivative with respect to W-.
            partialDeriv_[1] = findPartialMinus();

            // Solve the second semi-implicit equation for W- (j+1/2)
            WfieldsUpdate_[1] = stepWMinus(Wfields_[1], partialDeriv_[1]);;

            // Complete the full step by shifting the spatial average to zero
            shiftAverageZero(WfieldsUpdate_);

            // Update system and solve MDEs
            updateSystemFields();
            system().compute();
         }

      }

      // Failure: iteration counter itr reached maxItr without converging
      Log::file() << "Iterator failed to converge before the maximum number of iterations.\n";
      return 1;

   }

   template <int D>
   void SISIterator<D>::evaluateScatteringFnc()
   {
      // Evalute diblock scattering functions in Fourier space. This 
      // should only need to be done once as they are only functions of 
      // Fourier variable k, though this may change in a variable unit cell

      const IntVec<D> dftDim = gAA_.dftDimensions();
      FSArray<double,6> cellParam = system().unitCell().parameters();
      const double f = system().mixture().polymer(0).block(0).length() / 
                       system().mixture().polymer(0).length();
      
      // Determine actual dft mesh size, accounting for it being cut in half
      int dftSize = 1;
      for (int d = 0; d < D; d++) {
         dftSize*=dftDim[d];
      }
      
      for (int i = 0; i < dftSize; i++) {
         int temp = i;
         int sizeSlice = dftSize/dftDim[D-1];
         IntVec<D> coord;
         coord[D-1] = dftSize % sizeSlice;
         
         // For each dimension, determine coordinate in k-space array
         for (int d = D-2; d >= 0; d++) {
            temp -= coord[d+1]*sizeSlice;
            sizeSlice = dftSize/dftDim[d];
            coord[d] = temp % sizeSlice;
         }
         
         // find k norm by summing up over dimensions
         // Note: fftw.org says about the frequency at array position k 
         // that the actual frequnecy is equal to k/T where T is the 
         // periodicity in that dimension in real space. Because our 
         // arrays are periodic with the unit cell parameter in each
         // direction, I (Ryan) assumed that this was the cell parameter.

         double ksq = 0;
         for (int d = 0; d < D; d++) {
            ksq += pow((double)coord[d]/cellParam[d],2);
         }
         double k = pow(ksq,0.5);
         
         // Note: DFT fields are fields of fftw_complex numbers, with a real component [0] and imaginary component [1]
         gAA_[i][0] = 2/pow(k,4) * (f*pow(k,2) + exp(-pow(k,2)*f) - 1);
         gAB_[i][0] = 1/pow(k,4) * (1 - exp(-pow(k,2)*f))*(1 - exp(-pow(k,2)*(1-f)));
         gBB_[i][0] = 2/pow(k,4) * ( (1-f)*pow(k,2) + exp(-pow(k,2)*(1-f)) - 1);
      }
      
   }

   template <int D>
   void SISIterator<D>::shiftAverageZero(DArray<RField<D>> & fields)
   {
      // number of fields
      int nFields = fields.capacity();

      // for each field
      for (int i = 0; i < nFields; i++) {
         double avg = 0;
         int n = fields[i].capacity();
         
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
      // NOTE: Assumes that the MDEs were most recently solved on the 
      // system using the fields in Wfields_. Otherwise, this 
      // functional derivative is incorrect.

      // Get system data.
      const DArray<RField<D>> * CFields = &system().cFieldsRGrid();
      // Allocate temporary array
      RField<D> temp;
      temp.allocate((*CFields)[0].meshDimensions());
      // Compute functional derivative
      for (int i = 0; i < temp.capacity(); i++) { 
         temp[i] = (*CFields)[0][i] + (*CFields)[1][i] - 1;
      }

      return temp;
   }

   template <int D>
   RField<D> SISIterator<D>::findPartialMinus()
   {
      // NOTE: Assumes that the MDEs were most recently solved on the 
      // system using the fields in Wfields_. Otherwise, this 
      // functional derivative is incorrect.

      // Get system data
      const DArray<RField<D>> * CFields = &system().cFieldsRGrid();
      const double chiN = system().interaction().chi(0,1) * 
                           system().mixture().polymer(0).length();
      const double f = system().mixture().polymer(0).block(0).length() / 
                        system().mixture().polymer(0).length();
      // Allocate temporary array 
      RField<D> temp;
      temp.allocate((*CFields)[0].meshDimensions());
      // Compute functional derivative
      for (int i = 0; i < temp.capacity(); i++) { 
         temp[i] = (2*f - 1) + 2/chiN * Wfields_[1][i] - ((*CFields)[0][i] - (*CFields)[1][i]);
      }

      return temp;
   }

   template <int D>
   RField<D> SISIterator<D>::stepWPlus(const RField<D> & WPlus, const RField<D> & partialPlus)
   {
      // do FFT, solve in fourier space, FFT-inverse back
      // NOTE: How this (and its minus counterpart) are done will determine if
      // we need both Wfields_ and WfieldsUpdate_. Hopefully we don't!
      const IntVec<D> meshDim = system().mesh().dimensions();

      // FFT of current W+ field
      RFieldDft<D> WPlusDFT;
      WPlusDFT.allocate(meshDim);
      system().fft().forwardTransform(WPlus,WPlusDFT);

      // FFT of current functional derivative of H wrt W+
      RFieldDft<D> partialPlusDFT;
      partialPlusDFT.allocate(meshDim);
      system().fft().forwardTransform(partialPlus,partialPlusDFT);

      // FFT of updated field
      RFieldDft<D> WPlusUpdateDFT;
      WPlusUpdateDFT.allocate(meshDim);
      
      // Compute fourier transform of updated field
      // Determine actual dft mesh size, accounting for it being cut in half
      const IntVec<D> dftDim = WPlusDFT.dftDimensions();
      int dftSize = 1;
      for (int d = 0; d < D; d++) {
         dftSize*=dftDim[d];
      }
      // Note: DFT fields are fields of fftw_complex numbers, with a real component [0] and imaginary component [1]
      for (int i = 0; i < dftSize; i++) {
         WPlusUpdateDFT[i][0] = WPlusDFT[i][0] + dt_/(1 + dt_*(gAA_[i][0] + gBB_[i][0] + 2*gAB_[i][0] ))*partialPlusDFT[i][0];
      }

      RField<D> WPlusUpdate;
      WPlusUpdate.allocate(meshDim);
      system().fft().inverseTransform(WPlusUpdateDFT,WPlusUpdate);
      
      return WPlusUpdate;
   }

   template <int D>
   RField<D> SISIterator<D>::stepWMinus(const RField<D> & WMinus, const RField<D> & partialMinus)
   {
      // solve algebraically
      const IntVec<D> meshDim = system().mesh().dimensions();
      const double chiN = system().interaction().chi(0,1) * 
                           system().mixture().polymer(0).length();


      RField<D> WMinusUpdate;
      WMinusUpdate.allocate(meshDim);

      for (int i = 0; i < WMinusUpdate.capacity(); i++) {
         WMinusUpdate[i] = WMinus[i] - dt_/(1+dt_*2/chiN) * partialMinus[i];
      }

      return WMinusUpdate;
   }

   template <int D>
   DArray<RField<D>> SISIterator<D>::getWFields()
   {
      // get W fields, add/subtract to get W+ and W-, convert
      // to long vector format (can do we do this? should be able to)
      const int nx = Wfields_[0].capacity();
      const int nMon = Wfields_.capacity();

      DArray<RField<D>> Wfields;
      Wfields.allocate(nMon);
      for (int n = 0; n < nMon; n++) {
         Wfields[n].allocate(nx);
      }
      
      for (int i = 0; i < nx; i++) {
         Wfields[0][i] = 1/2*( system().wFieldRGrid(0)[i] + system().wFieldRGrid(1)[i] );
         Wfields[1][i] = 1/2*( system().wFieldRGrid(1)[i] - system().wFieldRGrid(0)[i] );
      }

      return Wfields;
      
   }

   template <int D>
   bool SISIterator<D>::isConverged() {

      double errorPlus, errorMinus, error;

      // Compute partial functional derivatives. Note: may be able to 
      // simplify this and do it fewer places! (they are already being computed
      // at each step of the algorithm...)
      RField<D> partialPlus = findPartialPlus();
      RField<D> partialMinus = findPartialMinus();

      // These should be zero if converged. Assess how close they are
      // to zero, depending on error type.

      // FOR NOW: default to maximum element. easy.
      errorPlus = 0;
      errorMinus = 0;
      for (int i = 0; i < partialPlus.capacity(); i++) {
         if (abs(partialPlus[i]) > errorPlus)
            errorPlus = abs(partialPlus[i]);
         if (abs(partialMinus[i]) > errorMinus) 
            errorMinus = abs(partialMinus[i]);

         if (errorMinus > errorPlus)
            error = errorMinus;
         else
            error = errorPlus;
      }
      
      if (errorMinus > errorPlus)
         error = errorMinus;
      else
         error = errorPlus;

      return error < epsilon_;
   }

   template <int D>
   void SISIterator<D>::updateSystemFields()
   {
      const int nx = Wfields_[0].capacity();
      const int nMon = Wfields_.capacity();

      DArray<typename System<D>::WField> Wupdate;

      Wupdate.allocate(nMon);
      for (int n = 0; n < nMon; n++) {
         Wupdate[n].allocate(nx);
      }
      for (int i = 0; i < nx; i++) {
         Wupdate[0][i] = WfieldsUpdate_[0][i] - WfieldsUpdate_[1][i];
         Wupdate[1][i] = WfieldsUpdate_[0][i] + WfieldsUpdate_[1][i];
      }
      system().setWRGrid(Wupdate);

   }
   

}
}
#endif
