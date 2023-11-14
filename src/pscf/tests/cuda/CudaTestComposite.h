#ifndef PSCF_TEST_CUDA_TEST_COMPOSITE_H
#define PSCF_TEST_CUDA_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

#include "CudaRandomTest.h"

TEST_COMPOSITE_BEGIN(CudaTestComposite)
TEST_COMPOSITE_ADD_UNIT(CudaRandomTest);
TEST_COMPOSITE_END

#endif
