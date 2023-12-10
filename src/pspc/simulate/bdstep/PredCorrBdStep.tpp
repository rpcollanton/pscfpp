#ifndef PSPC_PRED_CORR_BD_STEP_TPP
#define PSPC_PRED_CORR_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"

#include <pspc/simulate/BdSimulator.h>
#include <pspc/compressor/Compressor.h>
#include <pspc/System.h>
#include <pscf/math/IntVec.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   PredCorrBdStep<D>::PredCorrBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      wh_(),
      wf_(),
      eta_(),
      dwc_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   PredCorrBdStep<D>::~PredCorrBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void PredCorrBdStep<D>::readParameters(std::istream &in)
   {
      read<double>(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      wf_.allocate(nMonomer);
      wh_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         wf_[i].allocate(meshDimensions);
         wh_[i].allocate(meshDimensions);
      }
      eta_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         eta_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);

   }

   template <int D>
   void PredCorrBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(wf_.capacity() == nMonomer);
      UTIL_CHECK(wh_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(wf_[i].capacity() == meshSize);
         UTIL_CHECK(wh_[i].capacity() == meshSize);
      }
      UTIL_CHECK(eta_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(eta_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
   }

   template <int D>
   void PredCorrBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Copy current W fields from parent system
      for (i = 0; i < nMonomer; ++i) {
         wh_[i] = system().w().rgrid(i);
         wf_[i] = wh_[i];
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);

      // Construct all random displacements
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            eta[k] = b*random().gaussian();
         }
      }

      // Mid-step Predictor
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         RField<D> const & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            dwc_[k] = a*dc[k] + eta[k];
         }
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wn = wh_[i];
            evec = 0.5*simulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               wn[k] += evec*dwc_[k];
            }
         }
      }

      // Set modified fields at mid-point
      system().setWRGrid(wh_);
      system().compressor().compress();
      UTIL_CHECK(system().hasCFields());

      // Compute components and derivatives at mid-point
      simulator().clearData();
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();

      // Full step (corrector)
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         RField<D> const & eta = eta_[j];
         for (k = 0; k < meshSize; ++k) {
            dwc_[k] = a*dc[k] + eta[k];
         }
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wn = wf_[i];
            evec = simulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               wn[k] += evec*dwc_[k];
            }
         }
      }

      // Set fields at final point
      system().setWRGrid(wf_);
      system().compressor().compress();
      UTIL_CHECK(system().hasCFields());

      // Compute components and derivatives at final point
      simulator().clearData();
      simulator().computeWc();
      simulator().computeCc();
      simulator().computeDc();

   }

}
}
#endif
