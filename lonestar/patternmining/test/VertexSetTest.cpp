#include "galois/Galois.h"
#include "galois/Reduction.h"
#include "galois/Bag.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/ParallelSTL.h"

#include "galois/runtime/Profile.h"

#include "gtest/gtest.h"
#include "../VertexSet.h"


TEST(VertexSet, VertexSetInit) {
  uint32_t a[] = {0, 1, 2, 3, 4};
  uint32_t b[] = {2, 3};
  VertexSet<uint32_t,false,false> v0(4);
  VertexSet<uint32_t,false,false> v1(4);
  VertexSet<uint32_t,true,false>  v2(a,5,0);
  VertexSet<uint32_t,true,false>  v3(b,4,0);

  v2.differenceNoBound(v0, v3);
  v2.differenceWithBound(v1, v3, 3);

  uint32_t s0 = v2.differenceNoBoundSize(v3);
  uint32_t s1 = v2.differenceWithBoundSize(v3,3);

  v0.print();
  v1.print();
  std::cout << "s0: " << s0 << ", s1: " << s1 << std::endl;

  EXPECT_EQ (1, 1);
}
