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

namespace covariance_geometry
{

void Pose3DRPYTo3DQuaternion(const PoseRPY & pose_in, PoseQuaternion & pose_out)
{
  // Copy position
  pose_out.first = pose_in.first;
  // Axis-angle to quaternion
  RPYToQuaternion(pose_in.second, pose_out.second);
}

void Pose3DQuaternionTo3DRPY(const PoseQuaternion & pose_in, PoseRPY & pose_out)
{
  // Copy position
  pose_out.first = pose_in.first;
  // Quaternion to axis-angle
  QuaternionToRPY(pose_in.second, pose_out.second);
}

inline void RPYToQuaternion(const Eigen::Vector3d & rpy, Eigen::Quaterniond & quaternion)
{
  // Axis-angle to quaternion
  // Eigen::AngleAxisd rollAngle(rpy.x(), Eigen::Vector3d::UnitX());
  // Eigen::AngleAxisd pitchAngle(rpy.y(), Eigen::Vector3d::UnitY());
  // Eigen::AngleAxisd yawAngle(rpy.z(), Eigen::Vector3d::UnitZ());
  // quaternion = rollAngle * pitchAngle * yawAngle;

  // from:      https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  // quaternion.w() = cos(rpy.z() / 2) * cos(rpy.y() / 2) * cos(rpy.x() / 2) +
  //                  sin(rpy.z() / 2) * sin(rpy.y() / 2) * sin(rpy.x() / 2);
  // quaternion.x() = cos(rpy.z() / 2) * cos(rpy.y() / 2) * sin(rpy.x() / 2) -
  //                  sin(rpy.z() / 2) * sin(rpy.y() / 2) * cos(rpy.x() / 2);
  // quaternion.y() = cos(rpy.z() / 2) * sin(rpy.y() / 2) * cos(rpy.x() / 2) +
  //                  sin(rpy.z() / 2) * cos(rpy.y() / 2) * sin(rpy.x() / 2);
  // quaternion.z() = sin(rpy.z() / 2) * cos(rpy.y() / 2) * cos(rpy.x() / 2) -
  //                  cos(rpy.z() / 2) * sin(rpy.y() / 2) * sin(rpy.x() / 2);
  quaternion.w() = cos(rpy.x() / 2) * cos(rpy.y() / 2) * cos(rpy.z() / 2) +
    sin(rpy.x() / 2) * sin(rpy.y() / 2) * sin(rpy.z() / 2);
  quaternion.x() = sin(rpy.x() / 2) * cos(rpy.y() / 2) * cos(rpy.z() / 2) -
    cos(rpy.x() / 2) * sin(rpy.y() / 2) * sin(rpy.z() / 2);
  quaternion.y() = cos(rpy.x() / 2) * sin(rpy.y() / 2) * cos(rpy.z() / 2) +
    sin(rpy.x() / 2) * cos(rpy.y() / 2) * sin(rpy.z() / 2);
  quaternion.z() = cos(rpy.x() / 2) * cos(rpy.y() / 2) * sin(rpy.z() / 2) -
    sin(rpy.x() / 2) * sin(rpy.y() / 2) * cos(rpy.z() / 2);
}

inline void QuaternionToRPY(const Eigen::Quaterniond & quaternion, Eigen::Vector3d & rpy)
{
  // Quaternion to axis-angle
  // rpy = quaternion.toRotationMatrix().eulerAngles(0, 1, 2);

  // Equation 2.11 pag. 15 A tutorial on SE(3) transformation parameterizations and on-manifold optimization
  auto qx = quaternion.x();
  auto qy = quaternion.y();
  auto qz = quaternion.z();
  auto qw = quaternion.w();
  auto & roll = rpy.x();
  auto & pitch = rpy.y();
  auto & yaw = rpy.z();

  const auto discr = qw * qy - qx * qz;

  if (discr > 0.49999) {  // pitch = 90 deg
    pitch = 0.5 * M_PI;
    yaw = -2 * std::atan2(qx, qw);
    roll = 0;
  } else if (discr < -0.49999) {  // pitch =-90 deg
    pitch = -0.5 * M_PI;
    yaw = 2 * std::atan2(qx, qw);
    roll = 0;
  } else {  // Non-degenerate case:
    yaw = std::atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy * qy + qz * qz));
    pitch = std::asin(2 * discr);
    roll = std::atan2(2 * (qw * qx + qy * qz), 1 - 2 * (qx * qx + qy * qy));
  }
}
}  // namespace covariance_geometry
