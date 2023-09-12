// Copyright 2023 Andrea Ostuni, Giacomo Franchini - PIC4SeR
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "covariance_geometry/pose_composition.hpp"

#include <mrpt/poses/CPose3D.h>

#include "covariance_geometry/pose_representation.hpp"
#include "covariance_geometry/test_utils.hpp"
#include "gtest/gtest.h"

namespace covariance_geometry
{

TEST(PoseComposition, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1, identity, p_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  identity.first = {0, 0, 0};
  identity.second = {1, 0, 0, 0};  // W X Y Z
  ComposePose3DQuaternion(pq1, identity, p_out);
  EXPECT_TRUE(isApprox(pq1, p_out));
}

TEST(PoseComposition, HandleZeroPoseRPY)
{
  PoseRPY pr1, identity, p_out;
  pr1.first = {1.0, 0.0, 0.0};
  pr1.second = {0.9128709, 0.4082483, 0};
  identity.first = {0.0, 0.0, 0.0};
  identity.second = {0.0, 0.0, 0.0};
  ComposePose3DRPY(pr1, identity, p_out);
  EXPECT_TRUE(isApprox(pr1, p_out));
}

TEST(PoseComposition, HandleCompositionPoseQuaternion)
{
  PoseQuaternion pq1, pq2, pq3, p_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0.0, 0.0};
  pq2.first = {2, 0, 0};
  pq2.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  ComposePose3DQuaternion(pq1, pq2, p_out);

  mrpt::poses::CPose3D p1, p2;
  EigenToMRPT(pq1, p1);
  EigenToMRPT(pq2, p2);
  mrpt::poses::CPose3D p3 = p1 + p2;
  MRPTtoEigen(p3, pq3);

  EXPECT_TRUE(isApprox(pq3, p_out));
}

TEST(PoseComposition, HandleCompositionPoseRPY)
{
  PoseRPY pr1, pr2, pr3, p_out;
  pr1.first = {1, 0, 0};
  // pr1.second = {0.9128709, 0.4082483, 0};
  pr1.second = {1.5709, 0.0, 0.0};
  pr2.first = {2, 0, 0};
  // pr2.second = {0.9128709, 0.4082483, 0};
  pr2.second = {0.0, 1.5709, 0.0};
  ComposePose3DRPY(pr1, pr2, p_out);

  mrpt::poses::CPose3D p1, p2, p3;
  EigenToMRPT(pr1, p1);
  EigenToMRPT(pr2, p2);
  p3 = p1 + p2;
  MRPTtoEigen(p3, pr3);

  EXPECT_TRUE(isApprox(pr3, p_out));
}
}  // namespace covariance_geometry
