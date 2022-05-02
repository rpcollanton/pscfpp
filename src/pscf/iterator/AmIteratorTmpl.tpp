#ifndef PSCF_AM_ITERATOR_TMPL_TPP
#define PSCF_AM_ITERATOR_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/inter/ChiInteraction.h>
#include <util/containers/FArray.h>
#include <util/format/Dbl.h>
#include <util/misc/Timer.h>
#include <cmath>

namespace Pscf
{

   using namespace Util;

   /*
   * Constructor
   */
   template <typename Iterator, typename T>
   AmIteratorTmpl<Iterator,T>::AmIteratorTmpl()
    : epsilon_(0),
      lambda_(0),
      nHist_(0),
      maxHist_(0)
   {  setClassName("AmIteratorTmpl"); }

   /*
   * Destructor
   */
   template <typename Iterator, typename T>
   AmIteratorTmpl<Iterator,T>::~AmIteratorTmpl()
   {}

   /*
   * Read parameter file block.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::readParameters(std::istream& in)
   {
      errorType_ = "relNormResid"; // default type of error
      read(in, "maxItr", maxItr_);
      read(in, "epsilon", epsilon_);
      read(in, "maxHist", maxHist_);
      readOptional(in, "errorType", errorType_);

      if (!(errorType_ == "normResid" 
         || errorType_ == "maxResid"
         || errorType_ == "relNormResid")) {
         UTIL_THROW("Invalid iterator error type in parameter file.");
      }
      
   }

   /*
   * Setup and allocate memory required by iterator.
   */
   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::setup()
   {
      // Determine length of residual basis vectors
      nElem_ = nElements();

      // Allocate ring buffers
      fieldHists_.allocate(maxHist_+1);
      fieldBasis_.allocate(maxHist_+1);
      resHists_.allocate(maxHist_+1);
      resBasis_.allocate(maxHist_+1);
      
      // Allocate arrays used in iteration
      fieldTrial_.allocate(nElem_);
      resTrial_.allocate(nElem_);
      temp_.allocate(nElem_);
      
      // Allocate arrays/matrices used in coefficient calculation
      U_.allocate(maxHist_,maxHist_);
      v_.allocate(maxHist_);
      coeffs_.allocate(maxHist_);
   }

   /*
   * Solve iteratively.
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::solve()
   {

      // Preconditions:
      UTIL_CHECK(hasInitialGuess());

      // Timers for timing iteration components

      Timer timerMDE;
      Timer timerStress;
      Timer timerAM;
      Timer timerResid;
      Timer timerConverged;
      Timer timerCoeff;
      Timer timerOmega;
      Timer timerTotal;
      

      // Start overall timer 
      timerTotal.start();
      
      // Solve MDE for initial state
      timerMDE.start();
      evaluate();
      timerMDE.stop();

      // Iterative loop
      bool done;
      for (int itr = 0; itr < maxItr_; ++itr) {


         // Store current field in history ringbuffer
         getCurrent(temp_);
         fieldHists_.append(temp_);
         
         timerAM.start();

         Log::file()<<"---------------------"<<std::endl;
         Log::file()<<" Iteration  "<<itr<<std::endl;

         if (itr < maxHist_) {
            lambda_ = 1.0 - pow(0.9, itr+1);
            nHist_ = itr;
         } else {
            lambda_ = 1.0;
            nHist_ = maxHist_;
         }

         timerResid.start();
         computeResidual();
         timerResid.stop();

         // Test for convergence
         timerConverged.start();
         done = isConverged();
         timerConverged.stop();

         // Output details of this iteration to the log file
         outputToLog();

         if (done) {
            // Stop timers
            timerAM.stop();
            timerTotal.stop();

            Log::file() << "----------CONVERGED----------"<< std::endl;

            // Output timing results
            Log::file() << "\n";
            Log::file() << "Iterator times contributions:\n";
            Log::file() << "\n";
            Log::file() << "MDE solution:         " 
                        << timerMDE.time()  << " s,  "
                        << timerMDE.time()/timerTotal.time() << "\n";
            Log::file() << "residual computation: "  
                        << timerResid.time()  << " s,  "
                        << timerResid.time()/timerTotal.time() << "\n";
            Log::file() << "mixing coefficients:  "  
                        << timerCoeff.time()  << " s,  "
                        << timerCoeff.time()/timerTotal.time() << "\n";
            Log::file() << "checking convergence: "  
                        << timerConverged.time()  << " s,  "
                        << timerConverged.time()/timerTotal.time() << "\n";
            Log::file() << "updating guess:       "  
                        << timerOmega.time()  << " s,  "
                        << timerOmega.time()/timerTotal.time() << "\n";
            Log::file() << "total time:           "  
                        << timerTotal.time()   << " s  ";
            Log::file() << "\n\n";

            // Successful completion (i.e., converged within tolerance)
            cleanUp();
            return 0;

         } else {
            // Determine optimal linear combination coefficients for 
            // building the updated field guess
            timerCoeff.start();
            findResidCoeff();
            timerCoeff.stop();

            // Build the updated field
            timerOmega.start();
            updateGuess();
            timerOmega.stop();

            timerAM.stop();

            // Solve the fixed point equation
            timerMDE.start();
            evaluate();
            timerMDE.stop();

         }

      }
      // Failure: iteration counter itr reached maxItr without converging
      timerTotal.stop();
      
      Log::file() << "Iterator failed to converge before the maximum number of iterations.\n";
      cleanUp();
      return 1;
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::computeResidual()
   {
      // Get residuals
      getResidual(temp_);
      
      // Store residuals in residual history ringbuffer
      resHists_.append(temp_);

      return;
   }

   template <typename Iterator, typename T>
   bool AmIteratorTmpl<Iterator,T>::isConverged()
   {

      // Find max residual vector element
      double maxRes  = findMaxAbs(resHists_[0]);
      Log::file() << "Max Residual  = " << maxRes << std::endl;

      // Find norm of residual vector
      double normRes = findNorm(resHists_[0]);
      Log::file() << "Residual Norm = " << normRes << std::endl;

      // Find norm of residual vector relative to field
      double relNormRes = normRes/findNorm(fieldHists_[0]);
      Log::file() << "Relative Norm = " << relNormRes << std::endl;

      // Check if total error is below tolerance
      if (errorType_ == "normResid")
         return normRes < epsilon_;
      else if (errorType_ == "relNormResid")
         return relNormRes < epsilon_;
      else if (errorType_ == "maxResid")
         return maxRes < epsilon_;
      else
         UTIL_THROW("Invalid iterator error type in parameter file.");

   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::findResidCoeff()
   {
      // Initialize matrix and vector of residual dot products
      // if this is the first iteration
      if (nHist_ == 0) {
         for (int m = 0; m < maxHist_; ++m) {
            v_[m] = 0;
            coeffs_[m] = 0;
            for (int n = 0; n < maxHist_; ++n) {
               U_(m,n) = 0;
            }
         }
         return;
      }

      // Update basis vectors of residuals histories
      updateBasis(resBasis_, resHists_);

      // Update the U matrix and v vectors
      updateU(U_, resBasis_, nHist_);
      updateV(v_, resHists_[0], resBasis_, nHist_);

      // Solve matrix equation problem to get coefficients to minimize
      // the norm of the residual vector
      if (nHist_ == 1) {
         // Solve explicitly for coefficient
         coeffs_[0] = v_[0] / U_(0,0);
      } else
      if (nHist_ < maxHist_) { 
         // Create temporary smaller version of U_, v_, coeffs_
         // this is done to avoid reallocating U_ with each iteration.
         DMatrix<double> tempU;
         DArray<double> tempv,tempcoeffs;
         tempU.allocate(nHist_,nHist_);
         tempv.allocate(nHist_);
         tempcoeffs.allocate(nHist_);
         for (int i = 0; i < nHist_; ++i) {
            tempv[i] = v_[i];
            for (int j = 0; j < nHist_; ++j) {
               tempU(i,j) = U_(i,j);
            }
         }
         // Solve matrix equation
         LuSolver solver;
         solver.allocate(nHist_);
         solver.computeLU(tempU);
         solver.solve(tempv,tempcoeffs);
         // Transfer solution to full-sized member variable
         for (int i = 0; i < nHist_; ++i) {
            coeffs_[i] = tempcoeffs[i];
         }
      } else 
      if (nHist_ == maxHist_) {
         LuSolver solver;
         solver.allocate(maxHist_);
         solver.computeLU(U_);
         solver.solve(v_, coeffs_);
      }

      return;
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::updateGuess()
   {

      // Contribution of the last solution
      setEqual(fieldTrial_,fieldHists_[0]);
      setEqual(resTrial_, resHists_[0]);

      // If at least two histories
      if (nHist_ > 0) {
         // Update the basis vectors of field histories
         updateBasis(fieldBasis_,fieldHists_);
         // Combine histories into trial guess and predicted error
         addHistories(fieldTrial_, fieldBasis_, coeffs_, nHist_);
         addHistories(resTrial_, resBasis_, coeffs_, nHist_);
      }

      // Correct for predicted error
      addPredictedError(fieldTrial_,resTrial_,lambda_);

      // Send out updated guess to the iterator mediator
      update(fieldTrial_);
      
      return;
   }

   template <typename Iterator, typename T>
   void AmIteratorTmpl<Iterator,T>::cleanUp()
   {
      // Clear ring buffers
      resHists_.clear();
      resBasis_.clear();
      fieldHists_.clear();
      fieldBasis_.clear();
      return;
   }

}
#endif
