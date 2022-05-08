#ifndef PSPC_TEST_ITERATOR_TEST_COMPOSITE_H
#define PSPC_TEST_ITERATOR_TEST_COMPOSITE_H

#include <test/CompositeTestRunner.h>

// include the headers for individual tests
#include "SISIteratorTest.h"

TEST_COMPOSITE_BEGIN(IteratorTestComposite)
TEST_COMPOSITE_ADD_UNIT(SISIteratorTest)
TEST_COMPOSITE_END

#endif
