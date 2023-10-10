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
const int NUM_IT = 50;
const Eigen::Vector3d coord = {0.111, -1.24, 0.35};                               // x, y, z
const Eigen::Vector3d rpy = {0.9067432, 0.4055079, 0.1055943};                    // roll, pitch, yaw
const Eigen::Quaterniond quat = {0.8746791, 0.4379822, 0.1581314, 0.1345454};     // w, x, y, z

const Eigen::Vector3d rpy_gl = {0.12, M_PI_2, 0.34};                              // roll, pitch, yaw
const Eigen::Quaterniond quat_gl = {0.6884861, 0.1612045, 0.6884861, 0.1612045};  // w, x, y, z

TEST(PoseComposition, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1, identity, p_out;
  pq1.first = coord;
  pq1.second = quat;
  identity.first = Eigen::Vector3d::Zero();
  identity.second = Eigen::Quaterniond::Identity();
  ComposePose3DQuaternion(pq1, identity, p_out);
  EXPECT_TRUE(isApprox(pq1, p_out));
}

TEST(PoseComposition, HandleZeroPoseRPY)
{
  PoseRPY pr1, identity, p_out;
  pr1.first = coord;
  pr1.second = rpy;
  identity.first = Eigen::Vector3d::Zero();
  identity.second = Eigen::Vector3d::Zero();
  ComposePose3DRPY(pr1, identity, p_out);
  EXPECT_TRUE(isApprox(pr1, p_out));
}

TEST(PoseComposition, HandleCompositionPoseQuaternion)
{
  PoseQuaternion pq1, pq2, pq3, p_out;
  pq1.first = coord;
  pq1.second = quat;
  pq2.first = coord + Eigen::Vector3d::Random();
  pq2.second = Eigen::Quaterniond::UnitRandom();
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
  pr1.first = coord;
  pr1.second = rpy;
  pr2.first = coord + Eigen::Vector3d::Random();
  pr2.second = rpy + Eigen::Vector3d::Random();
  ComposePose3DRPY(pr1, pr2, p_out);

  mrpt::poses::CPose3D p1, p2, p3;
  EigenToMRPT(pr1, p1);
  EigenToMRPT(pr2, p2);
  p3 = p1 + p2;
  MRPTtoEigen(p3, pr3);

  EXPECT_TRUE(isApprox(pr3, p_out));
}

TEST(PoseComposition, HandleLoopCompositionPoseQuaternion)
{
  PoseQuaternion pq1, pq2, pq3, p_out;
  pq1.first = Eigen::Vector3d::Random();
  pq1.second = Eigen::Quaterniond::UnitRandom();
  pq2.first = Eigen::Vector3d::Random();
  pq2.second = Eigen::Quaterniond::UnitRandom();

  mrpt::poses::CPose3D p1, p2, p3;
  EigenToMRPT(pq1, p1);
  EigenToMRPT(pq2, p2);

  for (int i = 0; i < NUM_IT; i++) {
    ComposePose3DQuaternion(pq1, pq2, p_out);
    p3 = p1 + p2;
    MRPTtoEigen(p3, pq3);
    EXPECT_TRUE(isApprox(pq3, p_out));
    pq1 = pq2;
    pq2 = p_out;
    p1 = p2;
    p2 = p3;
  }
}
}  // namespace covariance_geometry
