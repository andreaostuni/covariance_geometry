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

#include <iostream>

namespace covariance_geometry
{

void Pose3DQuaternionCovarianceTo3DRPYCovariance(
  const PoseQuaternionCovariance & pose_quaternion_covariance,
  PoseRPYCovariance & pose_rpy_covariance)
{
  // Convert pose
  Pose3DQuaternionTo3DRPY(pose_quaternion_covariance.first, pose_rpy_covariance.first);
  // Convert covariance
  covariance3DQuaternionTo3DRPY(
    pose_quaternion_covariance.first.second, pose_quaternion_covariance.second,
    pose_rpy_covariance.second);
}

void Pose3DRPYCovarianceTo3DQuaternionCovariance(
  const PoseRPYCovariance & pose_rpy_covariance,
  PoseQuaternionCovariance & pose_quaternion_covariance)
{
  // Convert pose
  Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance.first);
  // Convert covariance
  covariance3DRPYTo3DQuaternion(
    pose_rpy_covariance.first.second, pose_rpy_covariance.second,
    pose_quaternion_covariance.second);
}

void Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(
  const PoseRPYCovariance & pose_rpy_covariance,
  PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy)
{
  // Convert pose
  Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance_rpy.first);
  // Copy covariance
  pose_quaternion_covariance_rpy.second = pose_rpy_covariance.second;
}

void Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(
  const PoseQuaternionCovariance & pose_quaternion_covariance,
  PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy)
{
  // Copy pose
  pose_quaternion_covariance_rpy.first = pose_quaternion_covariance.first;
  // Convert covariance
  covariance3DQuaternionTo3DRPY(
    pose_quaternion_covariance.first.second, pose_quaternion_covariance.second,
    pose_quaternion_covariance_rpy.second);
}

void Pose3DQuaternionCovarianceRPYTo3DRPYCovariance(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy,
  PoseRPYCovariance & pose_rpy_covariance)
{
  // Convert pose
  Pose3DQuaternionTo3DRPY(pose_quaternion_covariance_rpy.first, pose_rpy_covariance.first);
  // Copy covariance
  pose_rpy_covariance.second = pose_quaternion_covariance_rpy.second;
}

void Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy,
  PoseQuaternionCovariance & pose_quaternion_covariance)
{
  // Copy pose
  pose_quaternion_covariance.first = pose_quaternion_covariance_rpy.first;
  // Convert pose rotation from quat to RPY
  Eigen::Vector3d pose_rpy;
  QuaternionToRPY(pose_quaternion_covariance_rpy.first.second, pose_rpy);
  // Convert covariance
  covariance3DRPYTo3DQuaternion(
    pose_rpy, pose_quaternion_covariance_rpy.second, pose_quaternion_covariance.second);
}

PoseQuaternion InversePose(const PoseQuaternion & pose_in)
{
  PoseQuaternion pose_out;
  // Inverse of quaternion
  pose_out.second = pose_in.second.inverse();
  // Inverse of translation
  pose_out.first = pose_out.second * -pose_in.first;
  return pose_out;
}

PoseRPY InversePose(const PoseRPY & pose_in)
{
  // Easier to convert to quaternion, invert, and then convert back to RPY
  PoseQuaternion pose_q;
  PoseRPY pose_out;
  Pose3DRPYTo3DQuaternion(pose_in, pose_q);
  pose_q = InversePose(pose_q);
  Pose3DQuaternionTo3DRPY(pose_q, pose_out);
  return pose_out;
}

Eigen::Matrix7d inverseCovarianceQuaternion(
  const Eigen::Matrix7d & covariance_quaternion,
  const PoseQuaternion & pose)
{
  // Equation 6.3 pag. 34 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix7d cov_inv;
  Eigen::Matrix7d j_qi = Eigen::Matrix7d::Zero();

  // Equation 4.4 pag. 27 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  // j_qi top 3x7 block
  double dx = -pose.first.x();
  double dy = -pose.first.y();
  double dz = -pose.first.z();
  double qx = pose.second.x();
  double qy = pose.second.y();
  double qz = pose.second.z();
  double qw = pose.second.w();

  Eigen::Matrix<double, 3, 7> j_qi_top;
  Eigen::Matrix<double, 3, 4> j_qi_top_right;
  Eigen::Matrix4d j44;

  j_qi_top(0, 0) =  2.0 * (qy * qy + qz * qz) - 1.0;
  j_qi_top(0, 1) = -2.0 * (qw * qz + qx * qy);
  j_qi_top(0, 2) =  2.0 * (qw * qy - qx * qz);

  j_qi_top(1, 0) =  2.0 * (qw * qz - qx * qy);
  j_qi_top(1, 1) =  2.0 * (qx * qx + qz * qz) - 1.0;
  j_qi_top(1, 2) = -2.0 * (qw * qx + qy * qz);

  j_qi_top(2, 0) = -2.0 * (qw * qy + qx * qz);
  j_qi_top(2, 1) =  2.0 * (qw * qx - qy * qz);
  j_qi_top(2, 2) =  2.0 * (qx * qx + qy * qy) - 1.0;

  j_qi_top(0, 3) =  qy * dy + qz * dz;
  j_qi_top(0, 4) =  qx * dy - 2.0 * qy * dx - qw * dz;
  j_qi_top(0, 5) =  qx * dz + qw * dy - 2.0 * qz * dx;
  j_qi_top(0, 6) = -qy * dz + qz * dy;

  j_qi_top(1, 3) =  qy * dx - 2.0 * qx * dy + qw * dz;
  j_qi_top(1, 4) =  qx * dx + qz * dz;
  j_qi_top(1, 5) = -qw * dx - 2.0 * qz * dy + qy * dz;
  j_qi_top(1, 6) =  qx * dz - qz * dx;

  j_qi_top(2, 3) =  qz * dx - qw * dy - 2.0 * qx * dz;
  j_qi_top(2, 4) =  qz * dy + qw * dx - 2.0 * qy * dz;
  j_qi_top(2, 5) =  qx * dx + qy * dy;
  j_qi_top(2, 6) =  qy * dx - qx * dy;

  jacobianQuaternionNormalization(pose.second, j44);
  j_qi_top.block<3, 4>(0, 3) *= 2 * j44;
  j_qi.block<3, 7>(0, 0) = j_qi_top;

  // j_qi bottom right 4x4 block
  Eigen::DiagonalMatrix<double, 4> D(-1.0, -1.0, -1.0, 1.0);
  j_qi.block<4, 4>(3, 3) += D * j44;
  return j_qi * covariance_quaternion * j_qi.transpose();
}

Eigen::Matrix6d inverseCovarianceRPY(const Eigen::Matrix6d & covariance_rpy, const PoseRPY & pose)
{
  // Easier to convert covariance to quaternion, invert, and then convert back to RPY
  PoseQuaternion pose_quaternion;
  Eigen::Matrix7d covariance_quaternion;
  // Convert pose and covariance from RPY to quaternion
  Pose3DRPYTo3DQuaternion(pose, pose_quaternion);
  covariance3DRPYTo3DQuaternion(pose.second, covariance_rpy, covariance_quaternion);
  // Invert pose and covariance
  covariance_quaternion = inverseCovarianceQuaternion(covariance_quaternion, pose_quaternion);
  // TODO: pose inversion is redundant, it is already done in the inverseCovarianceQuaternion
  pose_quaternion = InversePose(pose_quaternion);
  // Convert back to RPY
  Eigen::Matrix6d covariance_rpy_out;
  covariance3DQuaternionTo3DRPY(pose_quaternion.second, covariance_quaternion, covariance_rpy_out);
  return covariance_rpy_out;
}

PoseQuaternionCovariance inversePose3DQuaternionCovarianceQuaternion(
  const PoseQuaternionCovariance & pose_quaternion_covariance)
{
  PoseQuaternionCovariance pose_quaternion_covariance_out;
  // Invert pose
  pose_quaternion_covariance_out.first = InversePose(pose_quaternion_covariance.first);
  // Invert covariance
  pose_quaternion_covariance_out.second = inverseCovarianceQuaternion(
    pose_quaternion_covariance.second, pose_quaternion_covariance.first);
  return pose_quaternion_covariance_out;
}

PoseRPYCovariance inversePose3DRPYCovarianceRPY(const PoseRPYCovariance & pose_rpy_covariance)
{
  PoseRPYCovariance pose_rpy_covariance_out;
  // Invert pose
  pose_rpy_covariance_out.first = InversePose(pose_rpy_covariance.first);
  // Invert covariance
  pose_rpy_covariance_out.second = inverseCovarianceRPY(
    pose_rpy_covariance.second, pose_rpy_covariance.first);
  return pose_rpy_covariance_out;
}

PoseQuaternionCovarianceRPY inversePose3DQuaternionCovarianceRPY(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy)
{
  PoseQuaternionCovarianceRPY pose_quaternion_covariance_rpy_out;
  // Invert pose
  pose_quaternion_covariance_rpy_out.first = InversePose(pose_quaternion_covariance_rpy.first);
  // Invert covariance
  // First convert pose in rpy, needed for the inverse of the covariance
  PoseRPY pose_rpy;
  Pose3DQuaternionTo3DRPY(pose_quaternion_covariance_rpy.first, pose_rpy);
  pose_quaternion_covariance_rpy_out.second = inverseCovarianceRPY(
    pose_quaternion_covariance_rpy.second, pose_rpy);
  return pose_quaternion_covariance_rpy_out;
}
}  // namespace covariance_geometry
