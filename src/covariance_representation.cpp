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

#include "covariance_geometry/covariance_representation.hpp"
// #include "covariance_geometry/pose_representation.hpp"

using PoseQuaternion = std::pair<Eigen::Vector3d, Eigen::Quaterniond>;
using PoseRPY = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
namespace covariance_geometry
{
void covariance3DRPYTo3DQuaternion(
  const Eigen::Vector3d & rpy, const Eigen::Matrix6d & covariance_rpy,
  Eigen::Matrix7d & covariance_quaternion)
{
  // Equation 2.8 pag. 14 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix7_6d jacobian = Eigen::Matrix7_6d::Zero();
  jacobian3DRPYTo3DQuaternion(rpy, jacobian);
  covariance_quaternion.noalias() = jacobian * covariance_rpy * jacobian.transpose();
}

void covariance3DQuaternionTo3DRPY(
  const Eigen::Quaterniond & quaternion, const Eigen::Matrix7d & covariance_quaternion,
  Eigen::Matrix6d & covariance_rpy)
{
  // Equation 2.12 pag. 16 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix6_7d jacobian = Eigen::Matrix6_7d::Zero();
  jacobian3DQuaternionTo3DRPY(quaternion, jacobian);
  covariance_rpy.noalias() = jacobian * covariance_quaternion * jacobian.transpose();
}

void jacobian3DQuaternionTo3DRPY(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix6_7d & jacobian)
{
  // Equation 2.13 pag. 16 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix3_4d j34_tmp;
  jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  jacobianQuaternionToRPY(quaternion, j34_tmp);
  jacobian.block<3, 4>(3, 3) = j34_tmp;
}

void jacobian3DRPYTo3DQuaternion(const Eigen::Vector3d & rpy, Eigen::Matrix7_6d & jacobian)
{
  // Equation 2.9a pag. 14 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  jacobian.block<3, 3>(0, 0).diagonal() = Eigen::Vector3d::Ones();
  jacobianRPYToQuaternion(rpy, jacobian.block<4, 3>(3, 3));
}

void jacobianQuaternionNormalization(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix4d & jacobian)
{
  // Equation 1.7 pag. 11 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();

  jacobian(0, 0) = qw * qw + qy * qy + qz * qz;
  jacobian(1, 1) = qw * qw + qx * qx + qz * qz;
  jacobian(2, 2) = qw * qw + qx * qx + qy * qy;
  jacobian(3, 3) = qx * qx + qy * qy + qz * qz;
  
  jacobian(0, 1) = -qx * qy;
  jacobian(0, 2) = -qx * qz;
  jacobian(0, 3) = -qx * qw;
  
  jacobian(1, 0) = -qx * qy;
  jacobian(1, 2) = -qy * qz;
  jacobian(1, 3) = -qy * qw;

  jacobian(2, 0) = -qx * qz;
  jacobian(2, 1) = -qy * qz;
  jacobian(2, 3) = -qz * qw;

  jacobian(3, 0) = -qx * qw;
  jacobian(3, 1) = -qy * qw;
  jacobian(3, 2) = -qz * qw;

  jacobian /= std::pow(quaternion.norm(), 3);
}

void jacobianRPYToQuaternion(const Eigen::Vector3d & rpy, Eigen::Ref<Eigen::Matrix4_3d> jacobian)
{
  // Equation 2.9b pag. 14 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  const double r2 = 0.5 * rpy.x();
  const double p2 = 0.5 * rpy.y();
  const double y2 = 0.5 * rpy.z();   
  const double ccc = cos(r2) * cos(p2) * cos(y2);
  const double ccs = cos(r2) * cos(p2) * sin(y2);
  const double csc = cos(r2) * sin(p2) * cos(y2);
  const double css = cos(r2) * sin(p2) * sin(y2);
  const double scc = sin(r2) * cos(p2) * cos(y2);
  const double scs = sin(r2) * cos(p2) * sin(y2);
  const double ssc = sin(r2) * sin(p2) * cos(y2);
  const double sss = sin(r2) * sin(p2) * sin(y2);

  // dqx()/d(rpy)
  jacobian(0, 0) =  (ccc + sss);
  jacobian(0, 1) = -(ssc + ccs);
  jacobian(0, 2) = -(csc + scs);

  // dqy()/d(rpy)
  jacobian(1, 0) = (ccs - ssc);
  jacobian(1, 1) = (ccc - sss);
  jacobian(1, 2) = (scc - css);

  // dqz()/d(rpy)
  jacobian(2, 0) = -(csc + scs);
  jacobian(2, 1) = -(css + scc);
  jacobian(2, 2) =  (ccc + sss);

  // dqw()/d(rpy)
  jacobian(3, 0) = (css - scc);
  jacobian(3, 1) = (scs - csc);
  jacobian(3, 2) = (ssc - ccs);

  jacobian *=  0.5;
}

void jacobianQuaternionToRPY(const Eigen::Quaterniond & quaternion, Eigen::Matrix3_4d & jacobian)
{
  // Equation 2.14 pag. 16 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  // d(rpy)()/d(quaternion) = d(rpy)()/d(quaternion_norm) * d(quaternion_norm)()/d(quaternion)
  // d(rpy)()/d(quaternion_norm) = jacobian_rpy_norm
  Eigen::Matrix3_4d jacobian_rpy_norm;
  jacobianQuaternionNormalToRPY(quaternion.normalized(), jacobian_rpy_norm);

  // d(quaternion_norm)()/d(quaternion) = jacobian_norm
  Eigen::Matrix4d jacobian_norm;
  jacobianQuaternionNormalization(quaternion, jacobian_norm);

  // d(rpy)()/d(quaternion) = d(rpy)()/d(quaternion_norm) * d(quaternion_norm)()/d(quaternion)
  jacobian.noalias() = jacobian_rpy_norm * jacobian_norm;
}

void jacobianQuaternionNormalToRPY(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix3_4d & jacobian)
{
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();
  const auto discr = qw * qy - qx * qz;

  if (discr > 0.49999) {
    // pitch = 90 deg
    jacobian = Eigen::Matrix3_4d::Zero();
    jacobian(2, 0) = -2 / (qw * ((qx * qx / qw * qw) + 1));
    jacobian(2, 3) = (2 * qx) / (qw * qw * ((qx * qx / qw * qw) + 1));
    return;
  } else if (discr < -0.49999) {
    // pitch = -90 deg
    jacobian = Eigen::Matrix3_4d::Zero();
    jacobian(2, 0) = 2 / (qw * ((qx * qx / qw * qw) + 1));
    jacobian(2, 3) = (-2 * qx) / (qw * qw * ((qx * qx / qw * qw) + 1));
    return;
  } else {
    // Non-degenerate case:
    jacobian(0, 0) =
      -((2 * qw) / (2 * qx * qx + 2 * qy * qy - 1) -
      (4 * qx * (2 * qw * qx + 2 * qy * qz)) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2)) /
      (std::pow((2 * qw * qx + 2 * qy * qz), 2) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2) + 1);
    jacobian(0, 1) =
      -((2 * qz) / (2 * qx * qx + 2 * qy * qy - 1) -
      (4 * qy * (2 * qw * qx + 2 * qy * qz)) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2)) /
      (std::pow((2 * qw * qx + 2 * qy * qz), 2) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2) + 1);
    jacobian(0, 2) =
      -(2 * qy) /
      ((std::pow((2 * qw * qx + 2 * qy * qz), 2) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2) +
      1) *
      (2 * qx * qx + 2 * qy * qy - 1));
    jacobian(0, 3) =
      -(2 * qx) /
      ((std::pow((2 * qw * qx + 2 * qy * qz), 2) / std::pow((2 * qx * qx + 2 * qy * qy - 1), 2) +
      1) *
      (2 * qx * qx + 2 * qy * qy - 1));

    jacobian(1, 0) = -(2 * qz) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 1) = (2 * qw) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 2) = -(2 * qx) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 3) = (2 * qy) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));

    jacobian(2, 0) =
      -(2 * qy) /
      ((std::pow((2 * qw * qz + 2 * qx * qy), 2) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2) +
      1) *
      (2 * qy * qy + 2 * qz * qz - 1));
    jacobian(2, 1) =
      -((2 * qx) / (2 * qy * qy + 2 * qz * qz - 1) -
      (4 * qy * (2 * qw * qz + 2 * qx * qy)) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2)) /
      (std::pow((2 * qw * qz + 2 * qx * qy), 2) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2) + 1);
    jacobian(2, 2) =
      -((2 * qw) / (2 * qy * qy + 2 * qz * qz - 1) -
      (4 * qz * (2 * qw * qz + 2 * qx * qy)) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2)) /
      (std::pow((2 * qw * qz + 2 * qx * qy), 2) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2) + 1);
    jacobian(2, 3) =
      -(2 * qz) /
      ((std::pow((2 * qw * qz + 2 * qx * qy), 2) / std::pow((2 * qy * qy + 2 * qz * qz - 1), 2) +
      1) *
      (2 * qy * qy + 2 * qz * qz - 1));

    return;
  }
}
}  // namespace covariance_geometry
