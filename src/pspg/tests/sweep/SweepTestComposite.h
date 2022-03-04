#ifndef PSPG_TEST_SWEEP_TEST_COMPOSITE_H
#define PSPG_TEST_SWEEP_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "RGridFieldStateTest.h"
#include "SweepTest.h"

TEST_COMPOSITE_BEGIN(SweepTestComposite)
TEST_COMPOSITE_ADD_UNIT(RGridFieldStateTest)
TEST_COMPOSITE_ADD_UNIT(SweepTest)
TEST_COMPOSITE_END

#endif
