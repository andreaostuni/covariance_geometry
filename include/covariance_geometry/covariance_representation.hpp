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

#ifndef COVARIANCE_GEOMETRY_COVARIANCE_REPRESENTATION_HPP
#define COVARIANCE_GEOMETRY_COVARIANCE_REPRESENTATION_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace Eigen
{
using Matrix7d = Eigen::Matrix<double, 7, 7>;
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Matrix6_7d = Eigen::Matrix<double, 6, 7>;
using Matrix7_6d = Eigen::Matrix<double, 7, 6>;
using Matrix3_4d = Eigen::Matrix<double, 3, 4>;
using Matrix4_3d = Eigen::Matrix<double, 4, 3>;
using Matrix3_7d = Eigen::Matrix<double, 3, 7>;
}  // namespace Eigen

using PoseQuaternion = std::pair<Eigen::Vector3d, Eigen::Quaterniond>;
using PoseRPY = std::pair<Eigen::Vector3d, Eigen::Vector3d>;

namespace covariance_geometry
{
/*
  / @brief Convert a covariance matrix from RPY representation to quaternion representation
  */
void covariance3DRPYTo3DQuaternion(
  const Eigen::Vector3d & rpy, const Eigen::Matrix6d & covariance_rpy,
  Eigen::Matrix7d & covariance_quaternion);
/*
  / @brief Convert a covariance matrix from quaternion representation to RPY representation
  */
void covariance3DQuaternionTo3DRPY(
  const Eigen::Quaterniond & quaternion, const Eigen::Matrix7d & covariance_quaternion,
  Eigen::Matrix6d & covariance_rpy);
/*
  / @brief jacobian of the transformation from 3D pose + quaternion to 3D pose + RPY
  */
void jacobian3DQuaternionTo3DRPY(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix6_7d & jacobian);

/*
  / @brief jacobian of the transformation from 3D pose + RPY to 3D pose + quaternion
  */
void jacobian3DRPYTo3DQuaternion(const Eigen::Vector3d & rpy, Eigen::Matrix7_6d & jacobian);

/*
  / @brief jacobian of the normalization of a quaternion
  */
void jacobianQuaternionNormalization(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix4d & jacobian);

/*
  / @brief jacobian of the transformation from quaternion to RPY
  */

void jacobianQuaternionToRPY(const Eigen::Quaterniond & quaternion, Eigen::Ref<Eigen::Matrix3_4d> jacobian);

/*
  / @brief jacobian of the transformation from RPY to quaternion
  */

void jacobianRPYToQuaternion(const Eigen::Vector3d & rpy, Eigen::Ref<Eigen::Matrix4_3d> jacobian);

/*
  / @brief jacobian of the transformation from normalized quaternion to RPY
  */

void jacobianQuaternionNormalToRPY(
  const Eigen::Quaterniond & quaternion, Eigen::Matrix3_4d & jacobian);

}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_COVARIANCE_REPRESENTATION_HPP
