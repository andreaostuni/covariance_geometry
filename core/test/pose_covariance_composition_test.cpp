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

#include "covariance_geometry/pose_covariance_composition.hpp"

#include <mrpt/poses/CPose3DPDFGaussian.h>
#include <mrpt/poses/CPose3DQuatPDFGaussian.h>

#include "covariance_geometry/test_utils.hpp"
#include "gtest/gtest.h"

using namespace mrpt::poses;
namespace covariance_geometry
{

const int NUM_IT = 50;
const Eigen::Vector3d coord = {10.0, 0.15465, -3.0};                               // x, y, z
const Eigen::Vector3d rpy = {0.9067432, 0.4055079, 0.1055943};                    // roll, pitch, yaw
const Eigen::Quaterniond quat = {0.8746791, 0.4379822, 0.1581314, 0.1345454};     // w, x, y, z

const Eigen::Vector3d rpy_gl = {0.12, M_PI_2, 0.34};                              // roll, pitch, yaw
const Eigen::Quaterniond quat_gl = {0.6884861, 0.1612045, 0.6884861, 0.1612045};  // w, x, y, z

TEST(PoseCovarianceComposition, HandleZeroPoseQuaternion)
{
  PoseQuaternionCovariance pq, identity, pq_out, pq_mrpt_out;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq.first.first = coord;
  pq.first.second = quat;
  pq.second = cov_eigen;
  identity.first.first = Eigen::Vector3d::Zero();
  identity.first.second = Eigen::Quaterniond::Identity();
  identity.second = Eigen::Matrix7d::Identity();
  ComposePoseQuaternionCovariance(pq, identity, pq_out);

  CPose3DQuatPDFGaussian p1, p2, p3;
  p1.mean = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  p1.cov = cov_mrpt;
  p2.mean = {0, 0, 0, {1, 0, 0, 0}};
  p2.cov = Eigen::Matrix7d::Identity();
  p3 = p1 + p2;
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(p3, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceComposition, HandleZeroPoseRPY)
{
  PoseRPYCovariance pr, identity, pr_out, pr_mrpt_out;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pr.first.first = coord;
  pr.first.second = rpy;
  pr.second = cov_eigen;
  identity.first.first = Eigen::Vector3d::Zero();
  identity.first.second = Eigen::Vector3d::Zero();
  identity.second = Eigen::Matrix6d::Identity();
  ComposePoseRPYCovarianceRPY(pr, identity, pr_out);

  CPose3DPDFGaussian p1, p2, p3;
  p1.mean = {coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x()};
  p1.cov = cov_mrpt;
  p2.mean = {0, 0, 0, 0, 0, 0};
  p2.cov = Eigen::Matrix6d::Identity();
  p3 = p1 + p2;
  MRPTtoEigen<CPose3DPDFGaussian, PoseRPYCovariance>(p3, pr_mrpt_out);
  permuteCovariance(pr_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseRPYCovariance>(pr_mrpt_out, pr_out));
}

TEST(PoseCovarianceComposition, HandleZeroPoseROS)
{
  PoseQuaternionCovarianceRPY pq, identity, pq_out, pq_mrpt_out;
  PoseRPY pr;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq.first.first = coord;
  pq.first.second = quat;
  pq.second = cov_eigen;
  identity.first.first = Eigen::Vector3d::Zero();
  identity.first.second = Eigen::Quaterniond::Identity();
  identity.second = Eigen::Matrix6d::Identity();
  ComposePoseQuaternionCovarianceRPY(pq, identity, pq_out);

  CPose3DQuat pq1, pq2, pq3;
  CPose3DPDFGaussian pr1, pr2, pr3;
  pq1 = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq2 = {0, 0, 0, {1, 0, 0, 0}};
  pq3 = pq1 + pq2;
  Pose3DQuaternionTo3DRPY(pq.first, pr);
  EigenToMRPT<PoseRPY, CPose3D>(pr, pr1.mean);
  pr1.cov = cov_mrpt;
  pr2.mean = {0, 0, 0, 0, 0, 0};
  pr2.cov = Eigen::Matrix6d::Identity();
  pr3 = pr1 + pr2;
  MRPTtoEigen<CPose3DPDFGaussian, PoseQuaternionCovarianceRPY>(pr3, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovarianceRPY>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceComposition, HandleCompositionPoseQuaternion)
{
  PoseQuaternionCovariance pq1, pq2, p_out, pq_mrpt_out;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq1.first.first = coord;
  pq1.first.second = quat;
  pq1.second = cov_eigen;
  pq2.first.first = coord;
  pq2.first.second = quat;
  pq2.second = cov_eigen;
  ComposePoseQuaternionCovariance(pq1, pq2, p_out);

  mrpt::poses::CPose3DQuatPDFGaussian p1, p2, p3;
  p1.mean = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  p1.cov = cov_mrpt;
  p2.mean = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  p2.cov = cov_mrpt;
  p3 = p1 + p2;
  MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(p3, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, p_out));
}

TEST(PoseCovarianceComposition, HandleCompositionPoseRPY)
{
  PoseRPYCovariance pr1, pr2, p_out, pr_mrpt_out;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pr1.first.first = coord;
  pr1.first.second = rpy;
  pr1.second = cov_eigen;

  pr2.first.first = coord;
  pr2.first.second = rpy;
  pr2.second = cov_eigen;

  ComposePoseRPYCovarianceRPY(pr1, pr2, p_out);

  mrpt::poses::CPose3DPDFGaussian p1, p2, p3;
  p1.mean = {coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x()};
  p1.cov = cov_mrpt;
  p2.mean = {coord.x(), coord.y(), coord.z(), rpy.z(), rpy.y(), rpy.x()};
  p2.cov = cov_mrpt;
  p3 = p1 + p2;
  MRPTtoEigen<CPose3DPDFGaussian, PoseRPYCovariance>(p3, pr_mrpt_out);
  permuteCovariance(pr_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseRPYCovariance>(pr_mrpt_out, p_out));
}

TEST(PoseCovarianceComposition, HandleCompositionPoseROS)
{
  PoseQuaternionCovarianceRPY p1, p2, identity, pq_out, pq_mrpt_out;
  PoseRPY pr;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6);
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  p1.first.first = coord;
  p1.first.second = quat;
  p1.second = cov_eigen;
  p2.first.first = coord;
  p2.first.second = quat;
  p2.second = cov_eigen;
  ComposePoseQuaternionCovarianceRPY(p1, p2, pq_out);

  CPose3DQuat pq1, pq2, pq3;
  CPose3DPDFGaussian pr1, pr2, pr3;
  pq1 = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq2 = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq3 = pq1 + pq2;
  Pose3DQuaternionTo3DRPY(p1.first, pr);
  EigenToMRPT<PoseRPY, CPose3D>(pr, pr1.mean);
  pr1.cov = cov_mrpt;
  pr2.mean = pr1.mean;
  pr2.cov = cov_mrpt;
  pr3 = pr1 + pr2;
  MRPTtoEigen<CPose3DPDFGaussian, PoseQuaternionCovarianceRPY>(pr3, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovarianceRPY>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceComposition, HandleCompositionPoseROSsmallcov)
{
  PoseQuaternionCovarianceRPY p1, p2, identity, pq_out, pq_mrpt_out;
  PoseRPY pr;
  Eigen::Matrix6d cov_mrpt = generateRandomCovariance(6) * 1e-41;
  Eigen::Matrix6d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  p1.first.first = coord;
  p1.first.second = quat;
  p1.second = cov_eigen;
  p2.first.first = coord;
  p2.first.second = quat;
  p2.second = cov_eigen * 1e26;
  ComposePoseQuaternionCovarianceRPY(p1, p2, pq_out);

  CPose3DQuat pq1, pq2, pq3;
  CPose3DPDFGaussian pr1, pr2, pr3;
  pq1 = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq2 = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  pq3 = pq1 + pq2;
  Pose3DQuaternionTo3DRPY(p1.first, pr);
  EigenToMRPT<PoseRPY, CPose3D>(pr, pr1.mean);
  pr1.cov = cov_mrpt;
  pr2.mean = pr1.mean;
  pr2.cov = cov_mrpt * 1e26;
  pr3 = pr1 + pr2;
  MRPTtoEigen<CPose3DPDFGaussian, PoseQuaternionCovarianceRPY>(pr3, pq_mrpt_out);
  permuteCovariance(pq_mrpt_out.second);

  EXPECT_TRUE(isApprox<PoseQuaternionCovarianceRPY>(pq_mrpt_out, pq_out));
}

TEST(PoseCovarianceComposition, HandleLoopCompositionPoseQuaternion)
{
  PoseQuaternionCovariance pq1, pq2, pq3, p_out, pq_mrpt_out;
  Eigen::Matrix7d cov_mrpt = generateRandomCovariance(7);
  Eigen::Matrix7d cov_eigen = cov_mrpt;
  permuteCovariance(cov_eigen);

  pq1.first.first = coord;
  pq1.first.second = quat;
  pq1.second = cov_eigen;
  pq2.first.first = coord;
  pq2.first.second = quat;
  pq2.second = cov_eigen;

  mrpt::poses::CPose3DQuatPDFGaussian p1, p2, p3;
  p1.mean = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  p1.cov = cov_mrpt;
  p2.mean = {coord.x(), coord.y(), coord.z(), {quat.w(), quat.x(), quat.y(), quat.z()}};
  p2.cov = cov_mrpt;

  for (int i = 0; i < NUM_IT; i++) {
    ComposePoseQuaternionCovariance(pq1, pq2, p_out);
    p3 = p1 + p2;
    MRPTtoEigen<CPose3DQuatPDFGaussian, PoseQuaternionCovariance>(p3, pq_mrpt_out);
    permuteCovariance(pq_mrpt_out.second);
    EXPECT_TRUE(isApprox<PoseQuaternionCovariance>(pq_mrpt_out, p_out));
    pq1 = pq2;
    pq2 = p_out;
    p1 = p2;
    p2 = p3;
  }
}
}  // namespace covariance_geometry
