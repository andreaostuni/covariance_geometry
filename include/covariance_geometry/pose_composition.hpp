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

#ifndef COVARIANCE_GEOMETRY_POSE_COMPOSITION_HPP
#define COVARIANCE_GEOMETRY_POSE_COMPOSITION_HPP

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <utility>

namespace covariance_geometry
{

using PoseQuaternion = std::pair<Eigen::Vector3d, Eigen::Quaterniond>;
using PoseRPY = std::pair<Eigen::Vector3d, Eigen::Vector3d>;

/*
  / @brief Pose composition for 6D poses in quaternion form
  */

void ComposePose3DQuaternion(const PoseQuaternion& a, const PoseQuaternion& b, PoseQuaternion& pose_out);

/*
  / @brief Pose composition for 6D poses in RPY form
  */

void ComposePose3DRPY(const PoseRPY& a, const PoseRPY& b, PoseRPY& pose_out);

}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_POSE_COMPOSITION_HPP
