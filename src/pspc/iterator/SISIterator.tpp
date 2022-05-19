#ifndef PSPC_SIS_ITERATOR_TPP
#define PSPC_SIS_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <math.h>                // for pow
#include <pscf/math/LuSolver.h>  // for finding inverse jacobian
#include "SISIterator.h"
#include <pspc/System.h>
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   template <int D>
   SISIterator<D>::SISIterator(System<D>& system)
   : Iterator<D>(system)
   {}

   template <int D>
   SISIterator<D>::~SISIterator()
   {}

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
      // flexible unit cell? 
      isFlexible_ = 0;
      scaleStress_ = 10;
      readOptional(in, "isFlexible", isFlexible_);
      readOptional(in, "scaleStress", scaleStress_);
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
      partialDeriv_.allocate(nField);
      for (int i = 0; i < nField; i++) {
         Wfields_[i].allocate(meshDim);
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

      // Do an initial stress relaxation with a high tolerance, 
      // if the unit cell is flexible.
      if (isFlexible_) {
         relaxUnitCell(epsilon_/scaleStress_*100);
      }

      // Iterative loop 
      bool fieldsConverged;
      for (int itr = 0; itr < maxItr_; ++itr) {
         getWFields(Wfields_);

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<" Iteration  "<<itr<<std::endl;
         
         // Compute error and test it with an isConverged function
         fieldsConverged = isConverged();

         // If the fields error is low enough
         if (fieldsConverged) {
            
            // Now that the fields are converged at for this unit cell, check the stress
            // to verify it is low enough and optimize cell parameters if cell is flexible.
            if (isFlexible_) {
               
               Log::file() << "Checking stress..." << std::endl;
               double stressTol = epsilon_/scaleStress_;
               bool stressConverged = isCellConverged(stressTol);
               
               if (!stressConverged) {
                  // Attempt to relax unit cell
                  relaxUnitCell(stressTol);
                  
                  // Continue field relaxations in the new cell by returning to for loop.
                  fieldsConverged = false;
                  Log::file()<< "Relaxing potential fields in new cell..." << std::endl;
                  continue;

               } else {
                  for (int i = 0; i < system().unitCell().nParameter(); i++) {
                     Log::file() << "Stress " << i <<  " = " << system().mixture().stress(i) << std::endl;
                  }
               }

            }
            
            Log::file() << "----------CONVERGED----------"<< std::endl;
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
            Wfields_[0] = stepWPlus(Wfields_[0], partialDeriv_[0]);

            // Update system with these updated fields and re-solve MDEs
            updateSystemFields(Wfields_);
            system().compute();

            // Find the functional derivative with respect to W-.
            partialDeriv_[1] = findPartialMinus(Wfields_[1]);

            // Solve the second semi-implicit equation for W- (j+1/2)
            Wfields_[1] = stepWMinus(Wfields_[1], partialDeriv_[1]);

            // Complete the full step by shifting the spatial average to zero
            shiftAverageZero(Wfields_);

            // Update system and solve MDEs
            updateSystemFields(Wfields_);
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
      
      // Expand out cell parameters because we need to use them in each direction
      FSArray<double,D> cellLengths;
      for (int i = 0; i < D; i++) {
         if (i >= cellParam.size()) {
            cellLengths.append(cellLengths[i-1]);
         } else {
            cellLengths.append(cellParam[i]);
         }
      }

      // Manually set k=0 components using limit of expresion as k goes to 0
      gAA_[0][0] = f*f;
      gAB_[0][0] = (1-f)*f;
      gBB_[0][0] = (1-f)*(1-f);
      
      // Determine actual dft mesh size, accounting for it being cut in half
      int dftSize = 1;
      for (int d = 0; d < D; d++) {
         dftSize*=dftDim[d];
      }

      for (int i = 1; i < dftSize; i++) {
         int temp = i;
         int sizeSlice = dftSize/dftDim[D-1];
         IntVec<D> coord;
         // last case is special
         if (sizeSlice != 1) 
            coord[D-1] = temp % sizeSlice;
         else 
            coord[D-1] = temp;

         // For each dimension, determine coordinate in k-space array
         for (int d = D-2; d >= 0; d--) {
            temp -= coord[d+1]*sizeSlice;
            sizeSlice = dftSize/dftDim[d];
            // last case is special 
            if (sizeSlice != 1) 
               coord[d] = temp % sizeSlice;
            else 
               coord[d] = temp;
         }
         
         // find k norm by summing up over dimensions
         // Note: fftw.org says about the frequency at array position k 
         // that the actual frequnecy is equal to k/T where T is the 
         // periodicity in that dimension in real space. Because our 
         // arrays are periodic with the unit cell parameter in each
         // direction, I (Ryan) assumed that this was the cell parameter.
         double ksq = 0;
         for (int d = 0; d < D; d++) {
            ksq += pow((double)coord[d]/cellLengths[d],2);
         }
         double k = pow(ksq,0.5);
         
         // Note: DFT fields are fields of fftw_complex numbers, with a real component [0] and imaginary component [1]
         // Divided by size to account for the factor of size that is picked up during the IFFT.
         gAA_[i][0] = 2/pow(k,4) * (f*pow(k,2) + exp(-pow(k,2)*f) - 1);
         gAA_[1][1] = 0.0;
         gAB_[i][0] = 1/pow(k,4) * (1 - exp(-pow(k,2)*f))*(1 - exp(-pow(k,2)*(1-f)));
         gAB_[1][1] = 0.0;
         gBB_[i][0] = 2/pow(k,4) * ( (1-f)*pow(k,2) + exp(-pow(k,2)*(1-f)) - 1);
         gBB_[1][1] = 0.0;
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
            fields[i][j] -= avg;
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
   RField<D> SISIterator<D>::findPartialMinus(const RField<D> & WMinus)
   {

      // Get system data
      const DArray<RField<D>> * CFields = &system().cFieldsRGrid();
      const double chiN = system().interaction().chi(0,1); // * system().mixture().polymer(0).length();
      const double f = system().mixture().polymer(0).block(0).length() / 
                        system().mixture().polymer(0).length();
      // Allocate temporary array 
      RField<D> temp;
      temp.allocate((*CFields)[0].meshDimensions());
      // Compute functional derivative
      for (int i = 0; i < temp.capacity(); i++) { 
         temp[i] = (2*f - 1) + 2/chiN * WMinus[i] - ((*CFields)[0][i] - (*CFields)[1][i]);
      }

      return temp;
   }

   template <int D>
   RField<D> SISIterator<D>::stepWPlus(const RField<D> & WPlus, const RField<D> & partialPlus)
   {
      // do FFT, solve in fourier space, FFT-inverse back

      const IntVec<D> meshDim = system().mesh().dimensions();
      const double size = WPlus.capacity();

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
         WPlusUpdateDFT[i][0] = WPlusDFT[i][0]*size + dt_/( 1 + dt_*(gAA_[i][0] + gBB_[i][0] + 2*gAB_[i][0]) )*partialPlusDFT[i][0]*size;
         WPlusUpdateDFT[i][0] /= size; // accounting for unnormalized fftw
         WPlusUpdateDFT[i][1] = WPlusDFT[i][1] + dt_*partialPlusDFT[i][1];
      }

      // Manually set 0th element of fourier transform to zero to set average to zero
      WPlusUpdateDFT[0][0] = 0.0; 

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
      const double chiN = system().interaction().chi(0,1); //* system().mixture().polymer(0).length();


      RField<D> WMinusUpdate;
      WMinusUpdate.allocate(meshDim);

      for (int i = 0; i < WMinusUpdate.capacity(); i++) {
         WMinusUpdate[i] = WMinus[i] - dt_/(1 + dt_*2/chiN) * partialMinus[i];
      }
      
      return WMinusUpdate;
   }

   template <int D>
   void SISIterator<D>::getWFields(DArray<RField<D>> & Wfields)
   {
      // get W fields, add/subtract to get W+ and W-, convert
      // to long vector format (can do we do this? should be able to)
      const int nx = Wfields[0].capacity();
      
      for (int i = 0; i < nx; i++) {
         Wfields[0][i] = 0.5*( system().wFieldRGrid(0)[i] + system().wFieldRGrid(1)[i] );
         Wfields[1][i] = 0.5*( system().wFieldRGrid(1)[i] - system().wFieldRGrid(0)[i] );
      }

      

      return;
   }

   template <int D> 
   DMatrix<double> SISIterator<D>::computeStressJacobian(FSArray<double,6> param)
   {
      const int nParam = system().unitCell().nParameter();

      // current and incremented stresses
      DArray<double> stresses, incrStresses;
      stresses.allocate(nParam);
      incrStresses.allocate(nParam);

      // get current stresses
      for (int i = 0; i < nParam; i++) {
         stresses[i] = system().mixture().stress(i);
      }

      // Parameter step size used for numerically approximating derivatives
      const double dparam = 1E-8;

      // Workspace for incrementing parameters
      FSArray<double,6> incrParam;

      // Initialize and allocate the Jacobian
      DMatrix<double> J;
      J.allocate(nParam,nParam);

      // Vary each parameter and compute the stress response, storing it in the jacobian
      for (int i = 0; i < nParam; i++) {
         // Reset parameters to values in currParam and increment one of them
         incrParam = param;
         incrParam[i] += dparam;
         // Update system with incremented parameter and compute stress
         system().setUnitCell(incrParam);
         system().compute();
         system().mixture().computeStress();
         // Get incremented stresses
         for (int j = 0; j < nParam; j++) {
            incrStresses[j] = system().mixture().stress(i);
         }
         // Numerically approximate derivative
         for (int j = 0; j < nParam; j++) {
            J(j,i) = 1/dparam * (incrStresses[j]-stresses[j]);
         }
      }

      // Reset system to original parameters!
      system().setUnitCell(param);
      system().compute();
      system().mixture().computeStress();

      return J;
   }

   template <int D>
   void SISIterator<D>::relaxUnitCell(double tol)
   {
      int nParam = system().unitCell().nParameter();
      DArray<double> stresses;
      stresses.allocate(nParam);

      // Compute the stress
      system().mixture().computeStress();
      
      // Output current stress values.
      Log::file()<<"---------------------"<<std::endl;
      Log::file() << "Relaxing Unit Cell. Initial Stress:" << std::endl;
      for (int i = 0; i < nParam; i++) {
         stresses[i] = system().mixture().stress(i);
         Log::file() << "Stress " << i <<  " = " << stresses[i] << std::endl;
      }
      Log::file() << std::endl;

      // If stress is not low enough.. iterate the unit cell parameters using Newton-Raphson (NR).
      bool stressConverged = isCellConverged(tol);
      if (!stressConverged) {

         // Get current unit cell parameters
         FSArray<double, 6> param = system().unitCell().parameters();
         // Allocate matrices for the jacobian and inverse jacobian 
         DMatrix<double> J, Jinv;
         J.allocate(nParam,nParam);
         Jinv.allocate(nParam,nParam);

         // Iterate to solve the stress equations for the optimal unit cell parameters         
         for (int stressItr = 0; stressItr < 50; stressItr++) {
            // Output iteration number
            Log::file()<<" Unit Cell Iteration "<<stressItr<<std::endl;

            // Compute jacobian
            J = computeStressJacobian(param);

            // Compute the inverse jacobian
            if (nParam == 1) {
               Jinv(0,0) = 1/J(0,0);
            } else {
               LuSolver solver;
               solver.allocate(nParam);
               solver.computeLU(J);
               solver.inverse(Jinv);
            }

            // Adjust unit cell parameters with one step of NR
            for (int i = 0; i < nParam; i++) {
               for (int j = 0; j < nParam; j++) {
                  param[i] -= Jinv(i,j)*stresses[j];
               }
            }
            // Update unit cell parameters on system
            system().setUnitCell(param);

            // Solve MDEs
            system().compute();

            // Recompute stress and check unit cell parameter convergence
            stressConverged = isCellConverged(tol);

            // Update stress and output stress/unit cell information
            for (int i = 0; i < nParam; i++) {
               stresses[i] = system().mixture().stress(i);
               Log::file()<<" -- Param  " << i << " = " << param[i] << std::endl;
               Log::file()<<" -- Stress " << i << " = " << stresses[i] << std::endl;
            }
            
            if (stressConverged) break;
         }

         // If the loop finishes and stress is still not converged, then the hardcoded max iters was reached.
         if (!stressConverged) {
            Log::file() << "Iterator reached maximum number of stress relaxation iterations." << std::endl;
         } else {
            Log::file() << "\nUnit cell relaxed." << std::endl;
         }         
         evaluateScatteringFnc();
      }

   }

   template <int D>
   bool SISIterator<D>::isCellConverged(double tol) 
   {
      bool stressConverged = true;
      system().mixture().computeStress();
      
      int nParam = system().unitCell().nParameter();
      DArray<double> stresses;
      stresses.allocate(nParam);

      for (int i = 0; i < nParam; i++) {
         stresses[i] = system().mixture().stress(i);
         if (abs(stresses[i]) > tol) { 
            stressConverged = false;
            break;
         }
      }

      return stressConverged;
   }

   template <int D>
   bool SISIterator<D>::isConverged() 
   {

      double errorPlus, errorMinus, error;

      // Compute partial functional derivatives. Note: may be able to 
      // simplify this and do it fewer places! (they are already being computed
      // at each step of the algorithm...)
      RField<D> partialPlus = findPartialPlus();
      RField<D> partialMinus = findPartialMinus(Wfields_[1]);

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

      Log::file() << "Max Error (Pressure) = " << errorPlus << std::endl;
      Log::file() << "Max Error (Exchange) = " << errorMinus << std::endl;
      
      if (errorMinus > errorPlus)
         error = errorMinus;
      else
         error = errorPlus;

      return error < epsilon_;
   }

   template <int D>
   void SISIterator<D>::updateSystemFields(const DArray<RField<D>> & WfieldsUpdate)
   {
      const int nx = Wfields_[0].capacity();
      const int nMon = Wfields_.capacity();

      DArray<typename System<D>::WField> WChain;
      WChain.allocate(nMon);

      for (int n = 0; n < nMon; n++) {
         WChain[n].allocate(nx);
      }
      for (int i = 0; i < nx; i++) {
         WChain[0][i] = WfieldsUpdate[0][i] - WfieldsUpdate[1][i];
         WChain[1][i] = WfieldsUpdate[0][i] + WfieldsUpdate[1][i];
      }
      system().setWRGrid(WChain);

   }
   
}
}
#endif
