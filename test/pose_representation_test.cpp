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

#include "covariance_geometry/pose_covariance_representation.hpp"
#include "covariance_geometry/test_utils.hpp"

#include <mrpt/poses/CPose3D.h>
#include <mrpt/poses/CPose3DQuat.h>

#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

namespace covariance_geometry
{
const int NUM_IT = 100;
const Eigen::Vector3d coord = {0.111, -1.24, 0.35};                             // x, y, z
const Eigen::Vector3d rpy = {0.9067432, 0.4055079, 0.1055943};                  // roll, pitch, yaw
const Eigen::Quaterniond quat = {0.8746791, 0.4379822, 0.1581314, 0.1345454};   // w, x, y, z

const Eigen::Vector3d rpy_gl = {0.0, M_PI_2, 0.0};                              // roll, pitch, yaw
const Eigen::Quaterniond quat_gl = {0.7070701, 0.0, 0.7071434, 0.0};            // w, x, y, z

TEST(EigenMRPT, EigenToMRPTRPY)
{
  PoseRPY pr1, pr_out;
  mrpt::poses::CPose3D mrpt_pose;
  pr1.first = coord;
  pr1.second = rpy;
  EigenToMRPT(pr1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pr_out);
  // GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pr1, pr_out));
}

TEST(EigenMRPT, MRPTtoEigenQuat)
{
  PoseQuaternion pq1, pq_out;
  mrpt::poses::CPose3D mrpt_pose;
  pq1.first = coord;
  pq1.second = quat;
  EigenToMRPT(pq1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pq_out);
  // GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox(pq1, pq_out));
}

TEST(PoseConversion, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1;
  PoseRPY pr1, p_out;
  pq1.first = Eigen::Vector3d::Zero();
  pq1.second = Eigen::Quaterniond::Identity();
  pr1.first = Eigen::Vector3d::Zero();
  pr1.second = Eigen::Vector3d::Zero();
  Pose3DQuaternionTo3DRPY(pq1, p_out);
  EXPECT_TRUE(isApprox(pr1, p_out));
}

TEST(PoseConversion, HandleZeroPoseRPY)
{
  PoseQuaternion pq1, p_out;
  PoseRPY pr1;
  pq1.first = Eigen::Vector3d::Zero();
  pq1.second = Eigen::Quaterniond::Identity();
  pr1.first = Eigen::Vector3d::Zero();
  pr1.second = Eigen::Vector3d::Zero();
  Pose3DRPYTo3DQuaternion(pr1, p_out);
  EXPECT_TRUE(isApprox(pq1, p_out));
}

TEST(PoseConversion, HandleConversionPoseQuaternion)
{
  PoseQuaternion pq;
  PoseRPY pr, pr_out;
  pq.first = coord;
  pq.second = quat;
  Pose3DQuaternionTo3DRPY(pq, pr_out);

  mrpt::poses::CPose3DQuat mrpt_pose_quat;
  mrpt::poses::CPose3D mrpt_pose_rpy;
  mrpt_pose_quat = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  mrpt_pose_rpy = mrpt::poses::CPose3D(mrpt_pose_quat);
  MRPTtoEigen<mrpt::poses::CPose3D, PoseRPY>(mrpt_pose_rpy, pr);

  EXPECT_TRUE(isApprox(pr, pr_out));
}

TEST(PoseConversion, HandleConversionPoseRPY)
{
  PoseQuaternion pq, pq_out;
  PoseRPY pr;
  pr.first = coord;
  pr.second = rpy;
  Pose3DRPYTo3DQuaternion(pr, pq_out);

  mrpt::poses::CPose3DQuat mrpt_pose_quat;
  mrpt::poses::CPose3D mrpt_pose_rpy;
  mrpt_pose_rpy = {coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x()};
  mrpt_pose_quat = mrpt::poses::CPose3DQuat(mrpt_pose_rpy);
  MRPTtoEigen<mrpt::poses::CPose3DQuat, PoseQuaternion>(mrpt_pose_quat, pq);

  EXPECT_TRUE(isApprox(pq, pq_out));
}

TEST(PoseConversion, HandleConversionPoseRPY_90Pitch)
{
  PoseQuaternion pq, pq_out;
  PoseRPY pr;
  pr.first = coord;
  pr.second = rpy_gl;
  Pose3DRPYTo3DQuaternion(pr, pq_out);

  mrpt::poses::CPose3DQuat mrpt_pose_quat;
  mrpt::poses::CPose3D mrpt_pose_rpy;
  mrpt_pose_rpy = {coord.x(), coord.y(), coord.z(), rpy_gl.z(), rpy_gl.y(), rpy_gl.x()};
  mrpt_pose_quat = mrpt::poses::CPose3DQuat(mrpt_pose_rpy);
  MRPTtoEigen<mrpt::poses::CPose3DQuat, PoseQuaternion>(mrpt_pose_quat, pq);

  EXPECT_TRUE(isApprox(pq, pq_out));

  pr.first = coord;
  pr.second = -rpy_gl;
  Pose3DRPYTo3DQuaternion(pr, pq_out);

  mrpt_pose_rpy = {coord.x(), coord.y(), coord.z(), -rpy_gl.z(), -rpy_gl.y(), -rpy_gl.x()};
  mrpt_pose_quat = mrpt::poses::CPose3DQuat(mrpt_pose_rpy);
  MRPTtoEigen<mrpt::poses::CPose3DQuat, PoseQuaternion>(mrpt_pose_quat, pq);

  EXPECT_TRUE(isApprox(pq, pq_out));
}

TEST(PoseConversion, HandleConversionPoseQuaternion_90Pitch)
{
  PoseQuaternion pq;
  PoseRPY pr, pr_out;
  pq.first = coord;
  pq.second = quat_gl;
  Pose3DQuaternionTo3DRPY(pq, pr_out);

  mrpt::poses::CPose3DQuat mrpt_pose_quat;
  mrpt::poses::CPose3D mrpt_pose_rpy;
  mrpt_pose_quat =
  {coord.x(), coord.y(), coord.z(), {quat_gl.w(), quat_gl.x(), quat_gl.y(), quat_gl.z()}};
  mrpt_pose_rpy = mrpt::poses::CPose3D(mrpt_pose_quat);
  MRPTtoEigen<mrpt::poses::CPose3D, PoseRPY>(mrpt_pose_rpy, pr);
  if (pr.second.isApprox(Eigen::Vector3d(M_PI, M_PI_2, M_PI), 1e-3)) {
    pr.second = Eigen::Vector3d(0.0, M_PI_2, 0.0);
  }
  EXPECT_TRUE(isApprox(pr, pr_out));

  pq.first = coord;
  pq.second = {quat_gl.w(), -quat_gl.x(), -quat_gl.y(), -quat_gl.z()};
  Pose3DQuaternionTo3DRPY(pq, pr_out);

  mrpt_pose_quat =
  {coord.x(), coord.y(), coord.z(), {quat_gl.w(), -quat_gl.x(), -quat_gl.y(), -quat_gl.z()}};
  mrpt_pose_rpy = mrpt::poses::CPose3D(mrpt_pose_quat);
  MRPTtoEigen<mrpt::poses::CPose3D, PoseRPY>(mrpt_pose_rpy, pr);
  if (pr.second.isApprox(Eigen::Vector3d(M_PI, -M_PI_2, M_PI), 1e-3)) {
    pr.second = Eigen::Vector3d(0.0, -M_PI_2, 0.0);
  }

  EXPECT_TRUE(isApprox(pr, pr_out));
}

TEST(PoseConversion, CyclicPoseConversionQuaternion)
{
  PoseQuaternion pq1, pq2;
  PoseRPY pr1;
  pq1.first = coord;
  pq1.second = quat;
  Pose3DQuaternionTo3DRPY(pq1, pr1);
  Pose3DRPYTo3DQuaternion(pr1, pq2);
  EXPECT_TRUE(isApprox(pq1, pq2));
}

TEST(PoseConversion, CyclicPoseConversionRPY)
{
  PoseQuaternion pq1;
  PoseRPY pr1, pr2;
  pr1.first = coord;
  pr1.second = rpy;
  Pose3DRPYTo3DQuaternion(pr1, pq1);
  Pose3DQuaternionTo3DRPY(pq1, pr2);
  EXPECT_TRUE(isApprox(pr1, pr2));
}

TEST(PoseInversion, InvertPoseQuaternion)
{
  PoseQuaternion pq, pq_inv, pq_out;
  mrpt::poses::CPose3DQuat pq_mrpt, pq_mrpt_inv;

  pq.first = coord;
  pq.second = quat;
  pq_inv = InversePose(pq);

  pq_mrpt = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq_mrpt_inv = -pq_mrpt;

  MRPTtoEigen<mrpt::poses::CPose3DQuat, PoseQuaternion>(pq_mrpt_inv, pq_out);
  EXPECT_TRUE(isApprox(pq_inv, pq_out));
}

TEST(PoseInversion, InvertPoseRPY)
{
  PoseRPY pr, pr_inv, pr_out;
  mrpt::poses::CPose3D pr_mrpt, pr_mrpt_inv;

  pr.first = coord;
  pr.second = rpy;
  pr_inv = InversePose(pr);

  pr_mrpt = {coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x()};
  pr_mrpt_inv = -pr_mrpt;

  MRPTtoEigen<mrpt::poses::CPose3D, PoseRPY>(pr_mrpt_inv, pr_out);
  EXPECT_TRUE(isApprox(pr_inv, pr_out));
}
}  // namespace covariance_geometry
