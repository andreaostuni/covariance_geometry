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

#include "covariance_geometry/pose_representation.hpp"

#include <mrpt/poses/CPose3D.h>

#include <cmath>
#include <iostream>

#include "covariance_geometry/test_utils.hpp"
#include "gtest/gtest.h"

namespace covariance_geometry
{

TEST(EigenMRPT, EigenToMRPTRPY)
{
  PoseRPY pr1, pr_out;
  mrpt::poses::CPose3D mrpt_pose;
  pr1.first = {1, 0, 0};
  pr1.second = {0.9128709, 0.4082483, 0};
  EigenToMRPT(pr1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pr_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pr1, pr_out));

  pr1.first = {1, 0, 0};
  pr1.second = {0.9067432, 0.4055079, 0.1055943};
  EigenToMRPT(pr1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pr_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pr1, pr_out));
}

TEST(EigenMRPT, MRPTtoEigenQuat)
{
  PoseQuaternion pq1, pq_out;
  mrpt::poses::CPose3D mrpt_pose;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  EigenToMRPT(pq1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pq_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  EigenToMRPT(pq1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pq_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pq1, pq_out));
}

TEST(PoseConversion, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1;
  PoseRPY pr1, p_out;
  pq1.first = {0, 0, 0};
  pq1.second = {1, 0, 0, 0};
  pr1.first = {0, 0, 0};
  pr1.second = {0, 0, 0};
  Pose3DQuaternionTo3DRPY(pq1, p_out);
  EXPECT_TRUE(isApprox(pr1, p_out));
}

TEST(PoseConversion, HandleZeroPoseRPY)
{
  PoseQuaternion pq1, p_out;
  PoseRPY pr1;
  pq1.first = {0, 0, 0};
  pq1.second = {1, 0, 0, 0};
  pr1.first = {0, 0, 0};
  pr1.second = {0, 0, 0};
  Pose3DRPYTo3DQuaternion(pr1, p_out);
  EXPECT_TRUE(isApprox(pq1, p_out));
}

TEST(PoseConversion, HandleConversionPoseQuaternion)
{
  PoseQuaternion pq1;
  PoseRPY pr1, pr_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8410687, 0.0, 0};
  Pose3DQuaternionTo3DRPY(pq1, pr_out);
  EXPECT_TRUE(isApprox(pr1, pr_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8545253, 0.1538007, 0.1742029};
  Pose3DQuaternionTo3DRPY(pq1, pr_out);
  EXPECT_TRUE(isApprox(pr1, pr_out));
}

TEST(PoseConversion, HandleConversionPoseRPY)
{
  PoseQuaternion pq1, pq_out;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8410687, 0.0, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8545253, 0.1538007, 0.1742029};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox(pq1, pq_out));
}

TEST(PoseConversion, HandleConversionPoseRPY_90Pitch)
{
  PoseQuaternion pq1, pq_out;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.7071068, 0.0, 0.7071068, 0.0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.0, M_PI_2, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.7071068, 0.0, -0.7071068, 0.0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.0, -M_PI_2, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox(pq1, pq_out));
}

TEST(PoseConversion, CyclicPoseConversionQuaternion)
{
  PoseQuaternion pq1, pq2;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  Pose3DQuaternionTo3DRPY(pq1, pr1);
  Pose3DRPYTo3DQuaternion(pr1, pq2);
  EXPECT_TRUE(isApprox(pq1, pq2));
}

TEST(PoseConversion, CyclicPoseConversionRPY)
{
  PoseQuaternion pq1;
  PoseRPY pr1, pr2;
  pr1.first = {1, 0, 0};
  pr1.second = {0.8545253, 0.1538007, 0.1742029};
  Pose3DRPYTo3DQuaternion(pr1, pq1);
  Pose3DQuaternionTo3DRPY(pq1, pr2);
  EXPECT_TRUE(isApprox(pr1, pr2));
}
}  // namespace covariance_geometry
