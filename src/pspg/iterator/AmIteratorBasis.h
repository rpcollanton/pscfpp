#ifndef PSPG_AM_ITERATOR_BASIS_H
#define PSPG_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspg
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspg implementation of the Anderson Mixing iterator.
   *
   * \ingroup Pspg_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis : public AmIteratorTmpl<Iterator<D>, DArray<double> >
   {

   public:

      /**
      * Constructor.
      *   
      * \param system parent system object
      */
      AmIteratorBasis(System<D>& system);

      /**
      * Destructor.
      */ 
      ~AmIteratorBasis();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      using AmIteratorTmpl<Iterator<D>, DArray<double> >::setup;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::solve;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::setClassName;
      using Iterator<D>::isFlexible;

   protected:

      using ParamComposite::readOptional;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;

   private:

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      /**
      * Find norm of a residual vector.
      */
      double findNorm(DArray<double> const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(DArray<double> const & hist);

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of pase residual or field vectors
      */
      void updateBasis(RingBuffer< DArray<double> > & basis, 
                       RingBuffer< DArray<double> > const & hists);

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<DArray<double> > const & resBasis, 
                             int m, int n);

      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent current residual vector
      * \param resBasis RingBuffer of residual basis vectors
      * \param m row index for element of the v vector
      */
      double computeVDotProd(DArray<double> const & resCurrent, 
                             RingBuffer<DArray<double> > const & resBasis, 
                             int m);

      /**
      * Update the U matrix.
      * 
      * \param U U matrix (dot products of residual basis vectors)
      * \param resBasis RingBuffer residual basis vectors
      * \param nHist current number of previous states
      */
      void updateU(DMatrix<double> & U, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void updateV(DArray<double> & v, 
                   DArray<double> const & resCurrent, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Set a vector equal to another. 
      * 
      * \param a the field to be set (LHS, result)
      * \param b the field for it to be set to (RHS, input)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Compute trial field so as to minimize L2 norm of residual.
      * 
      * \param trial resulting trial field (output)
      * \param basis RingBuffer of residual basis vectors.
      * \param coeffs coefficients of basis vectors
      * \param nHist number of prior states stored
      */
      void addHistories(DArray<double>& trial, 
                        RingBuffer<DArray<double> > const & basis, 
                        DArray<double> coeffs, int nHist);

      /**
      * Add predicted error to the trial field.
      * 
      * \param fieldTrial trial field (input/output)
      * \param resTrial predicted error for current trial field
      * \param lambda Anderson-Mixing mixing parameter 
      */
      void addPredictedError(DArray<double>& fieldTrial, 
                             DArray<double> const & resTrial, 
                             double lambda);

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
     
      /** 
      * Compute the number of elements in field or residual.
      */
      int nElements();

      /*
      * Get the current state of the system.
      *
      * \param curr current field vector (output)
      */
      void getCurrent(DArray<double>& curr);

      /**
      * Solve MDE for current state of system.
      */
      void evaluate();

      /**
      * Gets the residual vector from system.
      *  
      * \param curr current residual vector (output)
      */
      void getResidual(DArray<double>& resid);

      /**
      * Update the system with a new trial field vector.
      *
      * \param newGuess trial field configuration
      */
      void update(DArray<double>& newGuess);

      /**
      * Output relevant system details to the iteration log file.
      */
      void outputToLog();

      // --- Private member functions specific to this implementation --- 
      
      cudaReal findAverage(cudaReal * const field, int n);

   };

} // namespace Pspg
} // namespace Pscf
#endif