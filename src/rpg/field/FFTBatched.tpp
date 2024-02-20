#ifndef RPG_FFT_BATCHED_TPP
#define RPG_FFT_BATCHED_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FFTBatched.h"
#include <pscf/cuda/GpuResources.h>

namespace Pscf {
namespace Rpg
{

static __global__ void scaleComplexData(cudaComplex* data, cudaReal scale, int size) {
   
   //write code that will scale
   int nThreads = blockDim.x * gridDim.x;
   int startId = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startId; i < size; i += nThreads ) {
      data[i].x *= scale;
      data[i].y *= scale;
   }
   
}

static __global__ void scaleRealData(cudaReal* data, cudaReal scale, int size) {
   
   int nThreads = blockDim.x * gridDim.x;
   int startId = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startId; i < size; i += nThreads ) {
      data[i] *= scale;
   }
   
}

   using namespace Util;

   /*
   * Default constructor.
   */
   template <int D>
   FFTBatched<D>::FFTBatched()
    : meshDimensions_(0),
      rSize_(0),
      kSize_(0),
      fPlan_(0),
      iPlan_(0),
      isSetup_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   FFTBatched<D>::~FFTBatched()
   {
      if (fPlan_) {
         cufftDestroy(fPlan_);
      }
      if (iPlan_) {
         cufftDestroy(iPlan_);
      }
   }

   /*
   * Check and (if necessary) setup mesh dimensions.
   */
   template <int D>
   void FFTBatched<D>::setup(RField<D>& rField, RFieldDft<D>& kField)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);
      IntVec<D> rDimensions = rField.meshDimensions();
      IntVec<D> kDimensions = kField.meshDimensions();
      UTIL_CHECK(rDimensions == kDimensions);

      // Set Mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(rDimensions[i] > 0);
         meshDimensions_[i] = rDimensions[i];
         rSize_ *= rDimensions[i];
         if (i < D - 1) {
            kSize_ *= rDimensions[i];
         } else {
            kSize_ *= (rDimensions[i]/2 + 1);
         }
      }

      //turned off to allow batch solution of MDE
      //UTIL_CHECK(rField.capacity() == rSize_);
      //UTIL_CHECK(kField.capacity() == kSize_);


      // Make FFTW plans (explicit specializations)
      makePlans(rField, kField);

      isSetup_ = true;
   }

   template <int D>
   void FFTBatched<D>::setup(const IntVec<D>& rDim, const IntVec<D>& kDim, int batchSize)
   {
      // Preconditions
      UTIL_CHECK(!isSetup_);

      // Set Mesh dimensions
      rSize_ = 1;
      kSize_ = 1;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(rDim[i] > 0);
         meshDimensions_[i] = rDim[i];
         rSize_ *= rDim[i];
         if (i < D - 1) {
            kSize_ *= rDim[i];
         } else {
            kSize_ *= (rDim[i]/2 + 1);
         }
      }

      //turned off to allow batch solution of MDE
      //UTIL_CHECK(rField.capacity() == rSize_);
      //UTIL_CHECK(kField.capacity() == kSize_);


      // Make FFTW plans (explicit specializations)
      makePlans(rDim, kDim, batchSize);

      isSetup_ = true;
   }

   /*
   * Make plans for batch size of 2
   */
   template <int D>
   void FFTBatched<D>::makePlans(RField<D>& rField, RFieldDft<D>& kField)
   {
      int n[D];
      int rdist = 1;
      int kdist = 1;
      for(int i = 0; i < D; i++) {
         n[i] = rField.meshDimensions()[i];
         if( i == D - 1 ) {
            kdist *= n[i]/2 + 1;
         } else {
            kdist *= n[i];
         }
         rdist *= n[i];
         
      }

      #ifdef SINGLE_PRECISION
      if(cufftPlanMany(&fPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_R2C, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_C2R, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #else
      if(cufftPlanMany(&fPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_D2Z, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_Z2D, 2) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #endif      
   }

   /*
   * Make plans for variable batch size
   */

   template <int D>
   void FFTBatched<D>::makePlans(const IntVec<D>& rDim, const IntVec<D>& kDim, int batchSize)
   {
      int n[D];
      int rdist = 1;
      int kdist = 1;
      for(int i = 0; i < D; i++) {
         n[i] = rDim[i];
         if( i == D - 1 ) {
            kdist *= n[i]/2 + 1;
         } else {
            kdist *= n[i];
         }
         rdist *= n[i];
         
      }

      #ifdef SINGLE_PRECISION
      if(cufftPlanMany(&fPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_R2C, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, D, n, //plan, rank, n
                       NULL, 1, kdist, //inembed, istride, idist
                       NULL, 1, rdist, //onembed, ostride, odist
                       CUFFT_C2R, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #else
      if(cufftPlanMany(&fPlan_, D, n, //plan, rank, n
                       NULL, 1, rdist, //inembed, istride, idist
                       NULL, 1, kdist, //onembed, ostride, odist
                       CUFFT_D2Z, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      if(cufftPlanMany(&iPlan_, D, n, //plan, rank, n
                       NULL, 1, kdist, //inembed, istride, idist
                       NULL, 1, rdist, //onembed, ostride, odist
                       CUFFT_Z2D, batchSize) != CUFFT_SUCCESS) {
         std::cout<<"plan creation failed "<<std::endl;
         exit(1);
      }
      #endif
   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFTBatched<D>::forwardTransform(RField<D>& rField, RFieldDft<D>& kField)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(rSize_*2, nBlocks, nThreads);

      // Check dimensions or setup
      if (isSetup_) {
         UTIL_CHECK(rField.capacity() == 2*rSize_);
         UTIL_CHECK(kField.capacity() == 2*kSize_);
      } else {
         //UTIL_CHECK(0);
         //should never reach here in a parallel block. breaks instantly
         setup(rField, kField);
      }

      // Copy rescaled input data prior to work array
      cudaReal scale = 1.0/cudaReal(rSize_);
      //scale for every batch
      scaleRealData<<<nBlocks, nThreads>>>(rField.cField(), scale, rSize_ * 2);
      
      //perform fft
      #ifdef SINGLE_PRECISION
      if(cufftExecR2C(fPlan_, rField.cField(), kField.cField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #else
      if(cufftExecD2Z(fPlan_, rField.cField(), kField.cField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #endif

   }

   /*
   * Execute forward transform.
   */
   template <int D>
   void FFTBatched<D>::forwardTransform
   (cudaReal* rField, cudaComplex* kField, int batchSize)
   {
      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(rSize_*batchSize, nBlocks, nThreads);

      // Check dimensions or setup
      if (isSetup_) {
      } else {
         std::cerr<<"Need to setup before running transform"<<std::endl;
         exit(1);
      }

      // Copy rescaled input data prior to work array
      cudaReal scale = 1.0/cudaReal(rSize_);
      //scale for every batch
      scaleRealData<<<nBlocks, nThreads>>>(rField, scale, rSize_ * batchSize);
      
      //perform fft
      #ifdef SINGLE_PRECISION
      if(cufftExecR2C(fPlan_, rField, kField) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #else
      if(cufftExecD2Z(fPlan_, rField, kField) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: forward"<<std::endl;
         return;
      }
      #endif

   }

   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFTBatched<D>::inverseTransform(RFieldDft<D>& kField, RField<D>& rField)
   {
      if (!isSetup_) {
         //UTIL_CHECK(0);
         setup(rField, kField);
      }

      #ifdef SINGLE_PRECISION
      if(cufftExecC2R(iPlan_, kField.cField(), rField.cField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #else
      if(cufftExecZ2D(iPlan_, kField.cField(), rField.cField()) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #endif
   
   }


   /*
   * Execute inverse (complex-to-real) transform.
   */
   template <int D>
   void FFTBatched<D>::inverseTransform
   (cudaComplex* kField, cudaReal* rField, int batchSize)
   {
      if (!isSetup_) {
         std::cerr<<"Need to setup before running FFTbatchec"<<std::endl;
      }

      #ifdef SINGLE_PRECISION
      if(cufftExecC2R(iPlan_, kField, rField) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #else
      if(cufftExecZ2D(iPlan_, kField, rField) != CUFFT_SUCCESS) {
         std::cout<<"CUFFT error: inverse"<<std::endl;
         return;
      }
      #endif
   
   }

}
}
#endif
