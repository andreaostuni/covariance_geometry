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

namespace covariance_geometry
{
void covariance3DRPYTo3DQuaternion(
  const Eigen::Vector3d & rpy, const Eigen::Matrix6d & covariance_rpy,
  Eigen::Matrix7d & covariance_quaternion)
{
  // Equation 2.8 pag. 14 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix7_6d jacobian = Eigen::Matrix7_6d::Zero();
  jacobian3DRPYTo3DQuaternion(rpy, jacobian);
  covariance_quaternion = jacobian * covariance_rpy * jacobian.transpose();
}

void covariance3DQuaternionTo3DRPY(
  const Eigen::Quaterniond & quaternion, const Eigen::Matrix7d & covariance_quaternion,
  Eigen::Matrix6d & covariance_rpy)
{
  // Equation 2.12 pag. 16 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  Eigen::Matrix6_7d jacobian = Eigen::Matrix6_7d::Zero();
  jacobian3DQuaternionTo3DRPY(quaternion, jacobian);
  covariance_rpy = jacobian * covariance_quaternion * jacobian.transpose();
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
  Eigen::Matrix4_3d j43_tmp;
  jacobian.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  jacobianRPYToQuaternion(rpy, j43_tmp);
  jacobian.block<4, 3>(3, 3) = j43_tmp;
}

void jacobianQuaternionNormalization(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix4d & jacobian)
{
  // Equation 1.7 pag. 11 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  // Eigen::Quaterniond is quaternion in the form (x,y,z,w)
  jacobian = Eigen::Matrix4d::Zero();
  jacobian(0, 0) = quaternion.w() * quaternion.w() + quaternion.y() * quaternion.y() +
    quaternion.z() * quaternion.z();
  jacobian(1, 1) = quaternion.w() * quaternion.w() + quaternion.x() * quaternion.x() +
    quaternion.z() * quaternion.z();
  jacobian(2, 2) = quaternion.w() * quaternion.w() + quaternion.x() * quaternion.x() +
    quaternion.y() * quaternion.y();
  jacobian(3, 3) = quaternion.x() * quaternion.x() + quaternion.y() * quaternion.y() +
    quaternion.z() * quaternion.z();
  jacobian(0, 1) = -quaternion.x() * quaternion.y();
  jacobian(0, 2) = -quaternion.x() * quaternion.z();
  jacobian(0, 3) = -quaternion.x() * quaternion.w();

  jacobian(1, 0) = -quaternion.x() * quaternion.y();
  jacobian(1, 2) = -quaternion.y() * quaternion.z();
  jacobian(1, 3) = -quaternion.y() * quaternion.w();

  jacobian(2, 0) = -quaternion.x() * quaternion.z();
  jacobian(2, 1) = -quaternion.y() * quaternion.z();
  jacobian(2, 3) = -quaternion.z() * quaternion.w();

  jacobian(3, 0) = -quaternion.x() * quaternion.w();
  jacobian(3, 1) = -quaternion.y() * quaternion.w();
  jacobian(3, 2) = -quaternion.z() * quaternion.w();

  // auto norm_factor = 1.0 / std::pow(quaternion.norm(), 3);
  // jacobian = norm_factor * jacobian;
  jacobian = jacobian / std::pow(quaternion.norm(), 3);
}

void jacobianRPYToQuaternion(const Eigen::Vector3d & rpy, Eigen::Matrix4_3d & jacobian)
{
  // Equation 2.9b pag. 14 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  // rpy(0) = roll, rpy(1) = pitch, rpy(2) = yaw

  // const double cy = cos(rpy.z() / 2.0), sy = sin(rpy.z() / 2.0);
  // const double cp = cos(rpy.y() / 2.0), sp = sin(rpy.y() / 2.0);
  // const double cr = cos(rpy.x() / 2.0), sr = sin(rpy.x() / 2.0);

  // const double ccc = cr * cp * cy;
  // const double ccs = cr * cp * sy;
  // const double csc = cr * sp * cy;
  // const double css = cr * sp * sy;
  // const double scc = sr * cp * cy;
  // const double scs = sr * cp * sy;
  // const double ssc = sr * sp * cy;
  // const double sss = sr * sp * sy;

  const double ccc = cos(rpy.x() / 2.0) * cos(rpy.y() / 2.0) * cos(rpy.z() / 2.0);
  const double ccs = cos(rpy.x() / 2.0) * cos(rpy.y() / 2.0) * sin(rpy.z() / 2.0);
  const double csc = cos(rpy.x() / 2.0) * sin(rpy.y() / 2.0) * cos(rpy.z() / 2.0);
  const double css = cos(rpy.x() / 2.0) * sin(rpy.y() / 2.0) * sin(rpy.z() / 2.0);
  const double scc = sin(rpy.x() / 2.0) * cos(rpy.y() / 2.0) * cos(rpy.z() / 2.0);
  const double scs = sin(rpy.x() / 2.0) * cos(rpy.y() / 2.0) * sin(rpy.z() / 2.0);
  const double ssc = sin(rpy.x() / 2.0) * sin(rpy.y() / 2.0) * cos(rpy.z() / 2.0);
  const double sss = sin(rpy.x() / 2.0) * sin(rpy.y() / 2.0) * sin(rpy.z() / 2.0);

  // dqx()/d(rpy)
  // jacobian(0, 0) = -(csc + scs);
  // jacobian(0, 1) = -(ssc + ccs);
  // jacobian(0, 2) = (ccc + sss);
  // -0.5 * (csc + scs), 0.5 * (-ssc - ccs), 0.5 * (ccc + sss),

  jacobian(0, 0) = 0.5 * (ccc + sss);
  jacobian(0, 1) = 0.5 * -(ssc + ccs);
  jacobian(0, 2) = 0.5 * -(csc + scs);

  // dqy()/d(rpy)
  // jacobian(1, 0) = (scc - css);
  // jacobian(1, 1) = (ccc - sss);
  // jacobian(1, 2) = (ccs - scc);
  // 0.5 * (scc - css),	 0.5 * (ccc - sss),	 0.5 * (ccs - ssc),

  jacobian(1, 0) = 0.5 * (ccs - ssc);
  jacobian(1, 1) = 0.5 * (ccc - sss);
  jacobian(1, 2) = 0.5 * (scc - css);

  // dqz()/d(rpy)
  // jacobian(2, 0) = (ccc + sss);
  // jacobian(2, 1) = -(css + scc);
  // jacobian(2, 2) = -(csc + scs);
  // 0.5 * (ccc + sss),	 0.5 * (-css - scc), -0.5 * (csc + scs);

  jacobian(2, 0) = 0.5 * -(csc + scs);
  jacobian(2, 1) = 0.5 * -(css + scc);
  jacobian(2, 2) = 0.5 * (ccc + sss);

  // dqw()/d(rpy)
  // jacobian(3, 0) = (ssc - ccs);
  // jacobian(3, 1) = (scs - csc);
  // jacobian(3, 2) = (css - scc);

  // -0.5 * (ccs - ssc), 0.5 * (-csc + scs), -0.5 * (scc - css),

  jacobian(3, 0) = 0.5 * (css - scc);
  jacobian(3, 1) = 0.5 * (scs - csc);
  jacobian(3, 2) = 0.5 * (ssc - ccs);
  // jacobian = 0.5 * jacobian;
  jacobian = (1e-6 < jacobian.array().abs()).select(jacobian, 0.0);
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
  jacobian = jacobian_rpy_norm * jacobian_norm;
}

void jacobianQuaternionNormalToRPY(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix3_4d & jacobian)
{
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();
  const auto discr = qw * qy - qx * qz;

  jacobian = Eigen::Matrix3_4d::Zero();
  if (discr > 0.49999) {  // pitch = 90 deg
    jacobian(0, 0) = -2 * qw / (qx * qx + qw * qw);
    jacobian(0, 3) = +2 * qx / (qx * qx + qw * qw);
    return;
  } else if (discr < -0.49999) {  // pitch = -90 deg
    jacobian(0, 0) = +2 * qw / (qx * qx + qw * qw);
    jacobian(0, 3) = -2 * qx / (qx * qx + qw * qw);
    return;
  } else {  // Non-degenerate case:
    // row 1 --> becomes row3:
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

    // row 2:
    jacobian(1, 0) = -(2 * qz) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 1) = (2 * qw) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 2) = -(2 * qx) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));
    jacobian(1, 3) = (2 * qy) / std::sqrt(1 - std::pow((2 * qw * qy - 2 * qx * qz), 2));

    // row 3 --> becomes row1:
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
    return;
  }
}

}  // namespace covariance_geometry
