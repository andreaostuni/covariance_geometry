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
#include "covariance_geometry/pose_composition.hpp"
#include "covariance_geometry/pose_covariance_representation.hpp"

#include <iostream>
namespace covariance_geometry
{

void ComposePoseRPYCovarianceRPY(
  const PoseRPYCovariance & a, const PoseRPYCovariance & b, PoseRPYCovariance & out)
{
  // Make a way around them and consider instead this path:
  //
  //      X(6D)       U(6D)
  //        |           |
  //        v           v
  //      X(7D)       U(7D)
  //        |           |
  //        +--- (+) ---+
  //              |
  //              v
  //            RES(7D)
  //              |
  //              v
  //            RES(6D)

  // Poses conversion
  PoseQuaternionCovariance a_q, b_q, out_q;
  Pose3DRPYCovarianceTo3DQuaternionCovariance(a, a_q);
  Pose3DRPYCovarianceTo3DQuaternionCovariance(b, b_q);
  // Pose composition with covariance in quaternion
  ComposePoseQuaternionCovariance(a_q, b_q, out_q);
  // Pose conversion
  Pose3DQuaternionCovarianceTo3DRPYCovariance(out_q, out);
}

void ComposePoseQuaternionCovarianceRPY(
  const PoseQuaternionCovarianceRPY & a, const PoseQuaternionCovarianceRPY & b,
  PoseQuaternionCovarianceRPY & out)
{
  //  Make a way around them and consider instead this path:
  //
  // X(7D) COV(RPY) U(7D) COV(RPY)
  //        |           |
  //        v           v
  //      X(7D)       U(7D)
  //        |           |
  //        +--- (+) ---+
  //              |
  //              v
  //            RES(7D)
  //              |
  //              v
  //       RES(7D) COV(RPY)
  //

  PoseQuaternionCovariance a_q, b_q, out_q;
  // Poses conversion
  Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(a, a_q);
  Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(b, b_q);
  // Pose composition with covariance in quaternion
  ComposePoseQuaternionCovariance(a_q, b_q, out_q);
  // Pose conversion
  Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(out_q, out);
}

void ComposePoseQuaternionCovariance(
  const PoseQuaternionCovariance & a, const PoseQuaternionCovariance & b,
  PoseQuaternionCovariance & out)
{
  // Equation 5.7 pag. 31 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  // cov(a * b) = J_a * cov(a) * J_a^T + J_b * cov(b) * J_b^T
  auto cov_a = a.second;
  auto cov_b = b.second;
  auto cov_out = out.second;
  auto pose_a = a.first;
  auto pose_b = b.first;
  auto pose_out = out.first;

  // Pose composition
  ComposePose3DQuaternion(pose_a, pose_b, pose_out);

  // Covariance composition

  // derivative of the poses composition
  Eigen::Matrix7d j_fqc_a = Eigen::Matrix7d::Zero();
  Eigen::Matrix7d j_fqc_b = Eigen::Matrix7d::Zero();

  // Jacobian of the quaternion normalization function, computed in the output pose quaternion
  Eigen::Matrix4d jqn_out;
  jacobianQuaternionNormalization(pose_out.second, jqn_out);

  // d_fqc_ / d_a = [ d_fqr_ / d_a; 0 f(q1,q2)]
  JacobianPosePoseCompositionA(pose_a, pose_b, j_fqc_a);
  j_fqc_a.block<4, 4>(3, 3) = jqn_out * j_fqc_a.block<4, 4>(3, 3);

  // d_fqc_ / d_b = [ 0 d_fqr_ / d_b; 0 f(q1,q2)]
  JacobianPosePoseCompositionB(pose_a, j_fqc_b);
  j_fqc_b.block<4, 4>(3, 3) = jqn_out * j_fqc_b.block<4, 4>(3, 3);

  cov_out = j_fqc_a * cov_a * j_fqc_a.transpose() + j_fqc_b * cov_b * j_fqc_b.transpose();

  out.first = pose_out;
  out.second = cov_out;
}

void JacobianPosePoseCompositionA(
  const PoseQuaternion & pose_a, const PoseQuaternion & pose_b, Eigen::Matrix7d & jacobian)
{
  auto qx_b = pose_b.second.x();
  auto qy_b = pose_b.second.y();
  auto qz_b = pose_b.second.z();
  auto qw_b = pose_b.second.w();

  // Equation 5.8 pag. 31 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix3_7d j37_temp;

  JacobianPosePointComposition(pose_a, pose_b.first, j37_temp);
  jacobian.block<3, 7>(0, 0) = j37_temp;

  jacobian(3, 3) = qw_b;
  jacobian(3, 4) = qz_b;
  jacobian(3, 5) = -qy_b;
  jacobian(3, 6) = qx_b;

  jacobian(4, 3) = -qz_b;
  jacobian(4, 4) = qw_b;
  jacobian(4, 5) = qx_b;
  jacobian(4, 6) = qy_b;

  jacobian(5, 3) = qy_b;
  jacobian(5, 4) = -qx_b;
  jacobian(5, 5) = qw_b;
  jacobian(5, 6) = qz_b;

  jacobian(6, 3) = -qx_b;
  jacobian(6, 4) = -qy_b;
  jacobian(6, 5) = -qz_b;
  jacobian(6, 6) = qw_b;
}

void JacobianPosePoseCompositionB(const PoseQuaternion & pose_a, Eigen::Matrix7d & jacobian)
{
  auto qx_a = pose_a.second.x();
  auto qy_a = pose_a.second.y();
  auto qz_a = pose_a.second.z();
  auto qw_a = pose_a.second.w();

  // Equation 5.9 pag. 31 A tutorial on SE(3) transformation parameterizations and on-manifold optimization

  Eigen::Matrix3d j33_temp;

  JacobianPosePointComposition(pose_a, j33_temp);
  jacobian.block<3, 3>(0, 0) = j33_temp;

  jacobian(3, 3) = qw_a;
  jacobian(3, 4) = -qz_a;
  jacobian(3, 5) = qy_a;
  jacobian(3, 6) = qx_a;

  jacobian(4, 3) = qz_a;
  jacobian(4, 4) = qw_a;
  jacobian(4, 5) = -qx_a;
  jacobian(4, 6) = qy_a;

  jacobian(5, 3) = -qy_a;
  jacobian(5, 4) = qx_a;
  jacobian(5, 5) = qw_a;
  jacobian(5, 6) = qz_a;

  jacobian(6, 3) = -qx_a;
  jacobian(6, 4) = -qy_a;
  jacobian(6, 5) = -qz_a;
  jacobian(6, 6) = qw_a;
}

void JacobianPosePointComposition(
  const PoseQuaternion & pose, const Eigen::Vector3d & point, Eigen::Matrix3_7d & jacobian)
{
  // Equation 3.8 pag. 24 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix3_4d j34_temp;
  jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  JacobianQuaternionPointComposition(pose.second, point, j34_temp);
  jacobian.block<3, 4>(0, 3) = j34_temp;
}

void JacobianPosePointComposition(const PoseQuaternion & pose, Eigen::Matrix3d & jacobian)
{
  // Equation 3.10 pag. 24 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  auto qx = pose.second.x();
  auto qy = pose.second.y();
  auto qz = pose.second.z();
  auto qw = pose.second.w();

  jacobian(0, 0) = 0.5 - qy * qy - qz * qz;
  jacobian(0, 1) = qx * qy - qw * qz;
  jacobian(0, 2) = qx * qz + qw * qy;

  jacobian(1, 0) = qx * qy + qw * qz;
  jacobian(1, 1) = 0.5 - qx * qx - qz * qz;
  jacobian(1, 2) = qy * qz - qw * qx;

  jacobian(2, 0) = qx * qz - qw * qy;
  jacobian(2, 1) = qw * qx + qy * qz;
  jacobian(2, 2) = 0.5 - qx * qx - qy * qy;
  jacobian = 2 * jacobian;
}

void JacobianQuaternionPointComposition(
  const Eigen::Quaterniond & quaternion, const Eigen::Vector3d & point,
  Eigen::Matrix3_4d & jacobian)
{
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();

  auto ax = point.x();
  auto ay = point.y();
  auto az = point.z();

  // Equation 3.9 pag. 24 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  jacobian(0, 0) = qy * ay + qz * az;
  jacobian(0, 1) = -2 * qy * ax + qx * ay + qw * az;
  jacobian(0, 2) = -2 * qz * ax - qw * ay + qx * az;
  jacobian(0, 3) = -qz * ay + qy * az;

  jacobian(1, 0) = qy * ax - 2 * qx * ay - qw * az;
  jacobian(1, 1) = qx * ax + qz * az;
  jacobian(1, 2) = qw * ax - 2 * qz * ay + qy * az;
  jacobian(1, 3) = qz * ax - qx * az;

  jacobian(2, 0) = qz * ax + qw * ay - 2 * qx * az;
  jacobian(2, 1) = -qw * ax + qz * ay - 2 * qy * az;
  jacobian(2, 2) = qx * ax + qy * ay;
  jacobian(2, 3) = -qy * ax + qx * ay;

  Eigen::Matrix4d jqn;
  jacobianQuaternionNormalization(quaternion, jqn);
  jacobian = 2 * jacobian * jqn;
}
}  // namespace covariance_geometry
