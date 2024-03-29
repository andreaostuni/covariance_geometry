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

#include "covariance_geometry/pose_representation.hpp"

#include <cmath>
#include <iostream>
namespace covariance_geometry
{

void Pose3DRPYTo3DQuaternion(const PoseRPY & pose_in, PoseQuaternion & pose_out)
{
  // Copy position
  pose_out.first = pose_in.first;
  // Convert RPY to quaternion
  RPYToQuaternion(pose_in.second, pose_out.second);
}

void Pose3DQuaternionTo3DRPY(const PoseQuaternion & pose_in, PoseRPY & pose_out)
{
  // Copy position
  pose_out.first = pose_in.first;
  // Quaternion to RPY
  QuaternionToRPY(pose_in.second, pose_out.second);
}

inline void RPYToQuaternion(const Eigen::Vector3d & rpy, Eigen::Quaterniond & quaternion)
{
  const double half_roll = 0.5 * rpy.x();
  const double half_pitch = 0.5 * rpy.y();
  const double half_yaw = 0.5 * rpy.z();
  quaternion.w() = cos(half_roll) * cos(half_pitch) * cos(half_yaw) + sin(half_roll) * sin(
    half_pitch) * sin(half_yaw);
  quaternion.x() = sin(half_roll) * cos(half_pitch) * cos(half_yaw) - cos(half_roll) * sin(
    half_pitch) * sin(half_yaw);
  quaternion.y() = cos(half_roll) * sin(half_pitch) * cos(half_yaw) + sin(half_roll) * cos(
    half_pitch) * sin(half_yaw);
  quaternion.z() = cos(half_roll) * cos(half_pitch) * sin(half_yaw) - sin(half_roll) * sin(
    half_pitch) * cos(half_yaw);
  quaternion.normalize();
}

inline void QuaternionToRPY(const Eigen::Quaterniond & quaternion, Eigen::Vector3d & rpy)
{
  // Equation 2.11 pag. 15 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();
  auto & roll = rpy.x();
  auto & pitch = rpy.y();
  auto & yaw = rpy.z();

  const auto discr = qw * qy - qx * qz;
  if (discr > 0.49999) {
    // pitch = 90 deg
    pitch = 0.5 * M_PI;
    yaw = -2.0 * std::atan2(qx, qw);
    roll = 0.0;
  } else if (discr < -0.49999) {
    // pitch =-90 deg
    pitch = -0.5 * M_PI;
    yaw = 2.0 * std::atan2(qx, qw);
    roll = 0.0;
  } else {
    // Non-degenerate case:
    yaw = std::atan2(2.0 * (qw * qz + qx * qy), 1.0 - 2.0 * (qy * qy + qz * qz));
    pitch = std::asin(2.0 * discr);
    roll = std::atan2(2.0 * (qw * qx + qy * qz), 1.0 - 2.0 * (qx * qx + qy * qy));
  }
}
}  // namespace covariance_geometry
