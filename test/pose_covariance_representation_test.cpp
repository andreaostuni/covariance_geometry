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

#include <mrpt/poses/CPose3DPDFGaussian.h>
#include <mrpt/poses/CPose3DQuatPDFGaussian.h>

#include <cmath>
#include <iostream>

#include "covariance_geometry/covariance_representation.hpp"
#include "covariance_geometry/pose_representation.hpp"
#include "covariance_geometry/test_utils.hpp"
#include "gtest/gtest.h"

using namespace mrpt::poses;
namespace covariance_geometry
{

Eigen::Vector3d coord = {1, 0, 0};                                       // x, y, z
Eigen::Vector3d rpy = {0.9067432, 0.4055079, 0.1055943};                 // roll, pitch, yaw
Eigen::Quaterniond quat = {0.9067432, 0.4055079, 0.1055943, 0.0472232};  // w, x, y, z

TEST(CovarianceEigenMRPT, EigenToMRPTRPY)
{
  PoseRPYCovariance pr, pr_out;
  pr.first.first = coord;
  pr.first.second = rpy;
  pr.second = generateRandomCovariance(6);

  CPose3DPDFGaussian mrpt_pose;

  EigenToMRPT<PoseRPYCovariance, CPose3DPDFGaussian>(pr, mrpt_pose);
  MRPTtoEigen<CPose3DPDFGaussian, PoseRPYCovariance>(mrpt_pose, pr_out);
  EXPECT_TRUE(isApprox<PoseRPYCovariance>(pr, pr_out));
}

TEST(CovarianceEigenMRPT, MRPTtoEigenQuat)
{
  PoseQuaternionCovariance pq, pq_out;
  pq.first.first = coord;
  pq.first.second = quat;
  pq.second = generateRandomCovariance(7);

  CPose3DQuatPDFGaussian mrpt_pose;

  EigenToMRPT<PoseQuaternionCovariance, CPose3DQuatPDFGaussian>(pq, mrpt_pose);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose, pq_out);
  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq, pq_out));
}

TEST(PoseCovarianceConversion, HandleConversionPoseQuaternionCovariance)
{
  PoseQuaternionCovariance pq;
  PoseRPYCovariance pr_out, pr_mrpt_out;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq.first.first = coord;
  pq.first.second = quat;
  pq.second = cov_eigen;

  Pose3DQuaternionCovarianceTo3DRPYCovariance(pq, pr_out);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;
  mrpt_pose_q.mean =
    CPose3DQuat(coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()});
  mrpt_pose_q.cov = cov_mrpt;

  mrpt_pose_r = CPose3DPDFGaussian(mrpt_pose_q);
  MRPTtoEigen<CPose3DPDFGaussian, PoseRPYCovariance>(mrpt_pose_r, pr_mrpt_out);
  permuteCovariance(pr_mrpt_out.second);
  EXPECT_TRUE(isApprox<PoseRPYCovariance>(pr_mrpt_out, pr_out));
}

TEST(PoseCovarianceConversion, HandleConversionPoseRPYCovariance)
{
  PoseQuaternionCovariance pq_out, pq_mrpt_out;
  PoseRPYCovariance pr;
  pr.first.first = coord;
  pr.first.second = rpy;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);
  pr.second = cov_eigen;

  Pose3DRPYCovarianceTo3DQuaternionCovariance(pr, pq_out);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;
  mrpt_pose_r.mean = CPose3D(coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x());
  mrpt_pose_r.cov = cov_mrpt;
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);
  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceConversion, HandleConversionPoseROSCovariance)
{
  PoseRPYCovariance pr;
  PoseQuaternionCovarianceRPY pq_out;
  PoseQuaternionCovariance pq_mrpt_out;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pr.first.first = coord;
  pr.first.second = rpy;
  pr.second = cov_eigen;

  Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(pr, pq_out);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;
  mrpt_pose_r.mean = CPose3D(coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x());
  mrpt_pose_r.cov = cov_mrpt;
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);

  EXPECT_TRUE(isApprox(pq_mrpt_out.first, pq_out.first));
  // NB: here we are not checking the covariance since they're of different size
}

TEST(PoseCovarianceConversion, HandleConversionPoseRPYCovariance_90Pitch)
{
  PoseQuaternionCovariance pq_out, pq_mrpt_out;
  PoseRPYCovariance pr;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pr.first.first = coord;
  pr.first.second = {0.0, M_PI_2, 0.0};
  pr.second = cov_eigen;

  Pose3DRPYCovarianceTo3DQuaternionCovariance(pr, pq_out);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;

  mrpt_pose_r.mean = CPose3D(coord.x(), coord.y(), coord.z(), 0.0, M_PI_2, 0.0);
  mrpt_pose_r.cov = cov_mrpt;
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq_out));

  pr.first.first = coord;
  pr.first.second = {0.0, -M_PI_2, 0.0};
  pr.second = cov_eigen;
  Pose3DRPYCovarianceTo3DQuaternionCovariance(pr, pq_out);

  mrpt_pose_r.mean = CPose3D(coord.x(), coord.y(), coord.z(), 0.0, -M_PI_2, 0.0);
  mrpt_pose_r.cov = cov_mrpt;
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);
  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceConversion, CyclicPoseCovarianceConversionQuaternion)
{
  PoseQuaternionCovariance pq1, pq2, pq_mrpt_out;
  PoseRPYCovariance pr;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq1.first.first = coord;
  pq1.first.second = quat;
  pq1.second = cov_eigen;
  Pose3DQuaternionCovarianceTo3DRPYCovariance(pq1, pr);
  Pose3DRPYCovarianceTo3DQuaternionCovariance(pr, pq2);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;
  mrpt_pose_q.mean =
    CPose3DQuat(coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()});
  mrpt_pose_q.cov = cov_mrpt;
  mrpt_pose_r = CPose3DPDFGaussian(mrpt_pose_q);
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq2));
}

TEST(PoseCovarianceConversion, CyclicPoseCovarianceConversionRPY)
{
  PoseQuaternionCovariance pq;
  PoseRPYCovariance pr1, pr2;
  Eigen::Matrix6d cov_eigen = generateRandomCovariance(6);
  pr1.first.first = coord;
  pr1.first.second = rpy;
  pr1.second = cov_eigen;

  Pose3DRPYCovarianceTo3DQuaternionCovariance(pr1, pq);
  Pose3DQuaternionCovarianceTo3DRPYCovariance(pq, pr2);

  EXPECT_TRUE(isApprox<PoseRPYCovariance>(pr1, pr2));
}

TEST(PoseCovarianceConversion, CyclicPoseCovarianceConversionROS)
{
  PoseQuaternionCovariance pq1, pq2, pq_mrpt_out;
  PoseQuaternionCovarianceRPY pr;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);
  pq1.first.first = coord;
  pq1.first.second = quat;
  pq1.second = cov_eigen;
  Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(pq1, pr);
  Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(pr, pq2);

  CPose3DQuatPDFGaussian mrpt_pose_q;
  CPose3DPDFGaussian mrpt_pose_r;
  mrpt_pose_q.mean =
    CPose3DQuat(coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()});
  mrpt_pose_q.cov = cov_mrpt;
  mrpt_pose_r = CPose3DPDFGaussian(mrpt_pose_q);
  mrpt_pose_q = CPose3DQuatPDFGaussian(mrpt_pose_r);
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(mrpt_pose_q, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq2));
}
}  // namespace covariance_geometry
