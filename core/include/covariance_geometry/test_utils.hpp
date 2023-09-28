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

#ifndef COVARIANCE_GEOMETRY_TEST_UTILS_HPP
#define COVARIANCE_GEOMETRY_TEST_UTILS_HPP
#define GTEST_HAS_FILE_SYSTEM 1

#include "covariance_geometry/pose_covariance_representation.hpp"
#include "covariance_geometry/pose_representation.hpp"
#include "covariance_geometry/pose_composition.hpp"

#include <mrpt/poses/CPose3DPDFGaussian.h>
#include <mrpt/poses/CPose3DQuatPDFGaussian.h>

#include <iostream>
#include <iomanip>

#include "gtest/gtest.h"

#define GTEST_COUT std::cerr << "[          ] [ INFO ] " << std::setprecision(10)

namespace covariance_geometry
{

template<typename T>
bool isApprox(const T & a, const T & b, double epsilon = 1e-6)
{
  bool pose_equality = isApprox(a.first, b.first, epsilon);
  GTEST_COUT << "a.covariance: "
             << "\n"
             << a.second << std::endl;
  GTEST_COUT << "b.covariance: "
             << "\n"
    //  << std::setprecision(10)
             << b.second << std::endl;
  bool cov_equality = a.second.isApprox(b.second, epsilon);
  return pose_equality && cov_equality;
}

template<>
bool isApprox(const PoseQuaternion & a, const PoseQuaternion & b, double epsilon)
{
  GTEST_COUT << "a.first: " << a.first.transpose() << "\t";
  std::cout << "a.second: " << a.second << std::endl;
  GTEST_COUT << "b.first: " << b.first.transpose() << "\t";
  std::cout << "b.second: " << b.second << std::endl;
  Eigen::Quaterniond qa = a.second;
  ensurePositiveRealPart(qa);
  return a.first.isApprox(b.first, epsilon) && qa.isApprox(b.second, epsilon);
}

template<>
bool isApprox(const PoseRPY & a, const PoseRPY & b, double epsilon)
{
  GTEST_COUT << "a.first: " << a.first.transpose() << "\t";
  std::cout << "a.second: " << a.second.transpose() << std::endl;
  GTEST_COUT << "b.first: " << b.first.transpose() << "\t";
  std::cout << "b.second: " << b.second.transpose() << std::endl;
  return a.first.isApprox(b.first, epsilon) && a.second.isApprox(b.second, epsilon);
}

template<typename T, typename U>
bool MRPTtoEigen(const T & mrpt_pose, U & eigen_pose)
{
  eigen_pose.second = mrpt_pose.cov.asEigen();
  return MRPTtoEigen(mrpt_pose.getPoseMean(), eigen_pose.first);
}

template<>
bool MRPTtoEigen(const mrpt::poses::CPose3D & mrpt_pose, PoseRPY & eigen_pose)
{
  eigen_pose.first.x() = mrpt_pose.x();
  eigen_pose.first.y() = mrpt_pose.y();
  eigen_pose.first.z() = mrpt_pose.z();
  eigen_pose.second.x() = mrpt_pose.roll();
  eigen_pose.second.y() = mrpt_pose.pitch();
  eigen_pose.second.z() = mrpt_pose.yaw();
  return true;
}

template<>
bool MRPTtoEigen(const mrpt::poses::CPose3D & mrpt_pose, PoseQuaternion & eigen_pose)
{
  eigen_pose.first = mrpt_pose.m_coords.asEigen();
  mrpt::math::CQuaternionDouble q;
  mrpt_pose.getAsQuaternion(q);
  eigen_pose.second = Eigen::Quaterniond(q.r(), q.x(), q.y(), q.z());
  return true;
}

template<>
bool MRPTtoEigen(const mrpt::poses::CPose3DQuat & mrpt_pose, PoseQuaternion & eigen_pose)
{
  eigen_pose.first = mrpt_pose.m_coords.asEigen();
  mrpt::math::CQuaternionDouble q = mrpt_pose.quat();
  eigen_pose.second = Eigen::Quaterniond(q.r(), q.x(), q.y(), q.z());
  return true;
}

template<typename T, typename U>
bool EigenToMRPT(const T & eigen_pose, U & mrpt_pose)
{
  mrpt_pose.cov = eigen_pose.second;
  return EigenToMRPT(eigen_pose.first, mrpt_pose.mean);
}

template<>
bool EigenToMRPT(const PoseRPY & eigen_pose, mrpt::poses::CPose3D & mrpt_pose)
{
  mrpt_pose = mrpt_pose.FromXYZYawPitchRoll(
    eigen_pose.first(0), eigen_pose.first(1), eigen_pose.first(2), eigen_pose.second(2),
    eigen_pose.second(1), eigen_pose.second(0));
  return true;
}

template<>
bool EigenToMRPT(const PoseQuaternion & eigen_pose, mrpt::poses::CPose3D & mrpt_pose)
{
  mrpt::math::CQuaternionDouble q{
    eigen_pose.second.w(), eigen_pose.second.x(), eigen_pose.second.y(), eigen_pose.second.z()};
  mrpt_pose = mrpt_pose.FromQuaternionAndTranslation(
    q, eigen_pose.first.x(), eigen_pose.first.y(), eigen_pose.first.z());
  return true;
}

template<>
bool EigenToMRPT(const PoseQuaternion & eigen_pose, mrpt::poses::CPose3DQuat & mrpt_pose)
{
  mrpt::math::CQuaternionDouble q{
    eigen_pose.second.w(), eigen_pose.second.x(), eigen_pose.second.y(), eigen_pose.second.z()};
  mrpt_pose =
    mrpt::poses::CPose3DQuat(eigen_pose.first.x(), eigen_pose.first.y(), eigen_pose.first.z(), q);
  return true;
}

bool permuteCovariance(Eigen::Matrix7d & covariance)
{
  Eigen::PermutationMatrix<4> perm;
  perm.indices() = {1, 2, 3, 0};
  // process the upper right 3x4 block
  covariance.block<3, 4>(0, 3) = covariance.block<3, 4>(0, 3) * perm;
  // process the lower left 4x3 block
  covariance.block<4, 3>(3, 0) = perm.transpose() * covariance.block<4, 3>(3, 0);
  // process the lower right 4x4 block
  covariance.block<4, 4>(3, 3) = perm.transpose() * covariance.block<4, 4>(3, 3) * perm;
  return true;
}

bool permuteCovariance(Eigen::Matrix6d & covariance)
{
  Eigen::PermutationMatrix<3> perm;
  perm.indices() = {2, 1, 0};
  // process the upper right 3x3 block
  covariance.block<3, 3>(0, 3) = covariance.block<3, 3>(0, 3) * perm;
  // process the lower left 3x3 block
  covariance.block<3, 3>(3, 0) = perm.transpose() * covariance.block<3, 3>(3, 0);
  // process the lower right 3x3 block
  covariance.block<3, 3>(3, 3) = perm.transpose() * covariance.block<3, 3>(3, 3) * perm;
  return true;
}

Eigen::MatrixXd generateRandomCovariance(const int size)
{
  Eigen::MatrixXd cov = Eigen::MatrixXd::Random(size, size);
  return cov.selfadjointView<Eigen::Upper>();
}
}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_TEST_UTILS_HPP
