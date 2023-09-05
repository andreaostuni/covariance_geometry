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

#ifndef COVARIANCE_GEOMETRY_POSE_REPRESENTATION_HPP
#define COVARIANCE_GEOMETRY_POSE_REPRESENTATION_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace covariance_geometry
{

using PoseQuaternion = std::pair<Eigen::Vector3d, Eigen::Quaterniond>;
using PoseRPY = std::pair<Eigen::Vector3d, Eigen::Vector3d>;

/*
  / @brief Convert a pose from RPY representation to quaternion representation
  */
void Pose3DRPYTo3DQuaternion(const PoseRPY& pose_in, PoseQuaternion& pose_out);

/*
  / @brief Convert a pose from quaternion representation to RPY representation
  */
void Pose3DQuaternionTo3DRPY(const PoseQuaternion& pose_in, PoseRPY& pose_out);

/*
  / @brief Convert a RPY representation to quaternion representation
  */
void RPYToQuaternion(const Eigen::Vector3d& rpy, Eigen::Quaterniond& quaternion);

/*
  / @brief Convert a quaternion representation to RPY representation
  */
void QuaternionToRPY(const Eigen::Quaterniond& quaternion, Eigen::Vector3d& rpy);

}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_POSE_REPRESENTATION_HPP
