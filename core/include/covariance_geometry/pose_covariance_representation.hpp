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

#ifndef COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP
#define COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP

#include "covariance_geometry/covariance_representation.hpp"
#include "covariance_geometry/pose_representation.hpp"

namespace covariance_geometry
{
using PoseQuaternionCovariance = std::pair<PoseQuaternion, Eigen::Matrix7d>;
using PoseRPYCovariance = std::pair<PoseRPY, Eigen::Matrix6d>;
using PoseQuaternionCovarianceRPY = std::pair<PoseQuaternion, Eigen::Matrix6d>;  // Ros convention

/*
  / @brief Convert a Pose with covariance matrix from quaternion representation to RPY representation
  */
void Pose3DQuaternionCovarianceTo3DRPYCovariance(
  const PoseQuaternionCovariance & pose_quaternion_covariance,
  PoseRPYCovariance & pose_rpy_covariance);

/*
  / @brief Convert a Pose with covariance matrix from RPY representation to quaternion representation
  */
void Pose3DRPYCovarianceTo3DQuaternionCovariance(
  const PoseRPYCovariance & pose_rpy_covariance,
  PoseQuaternionCovariance & pose_quaternion_covariance);

/*
  / @brief Convert a Pose with covariance matrix from RPY representation to quaternion representation but covariance
  still in RPY (ROS convention)
  */
void Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(
  const PoseRPYCovariance & pose_rpy_covariance,
  PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy);

/*
  / @brief Convert a Pose with covariance matrix from quaternion representation to quaternion representation but
  covariance still in RPY (ROS convention)
  */
void Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(
  const PoseQuaternionCovariance & pose_quaternion_covariance,
  PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy);

/*
  / @brief Convert a Pose with covariance matrix from quaternion representation with covariance in RPY (ROS
  representation) to RPY representation
  */
void Pose3DQuaternionCovarianceRPYTo3DRPYCovariance(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy,
  PoseRPYCovariance & pose_rpy_covariance);

/*
  / @brief Convert a Pose with covariance matrix from quaternion representation with covariance in RPY (ROS
  representation) to quaternion representation
  */
void Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy,
  PoseQuaternionCovariance & pose_quaternion_covariance);

/*
  / @brief Inverse of a pose in quaternion representation
 */
PoseQuaternion InversePose(const PoseQuaternion & pose_in);

/*
  / @brief Inverse of a pose in RPY representation
 */
PoseRPY InversePose(const PoseRPY & pose_in);

/*
  / @brief Inverse of a covariance matrix in quaternion representation
  */

Eigen::Matrix7d inverseCovarianceQuaternion(
  const Eigen::Matrix7d & covariance_quaternion,
  const PoseQuaternion & pose);

/*
  / @brief Inverse of a covariance matrix in RPY representation
  This accepts pose in quaternion representation, since it is assumed to be used together with inversePose
  */

Eigen::Matrix6d inverseCovarianceRPY(const Eigen::Matrix6d & covariance_rpy, const PoseRPY & pose);

/*
  / @brief Inverse of a pose with covariance in quaternion representation
 */
PoseQuaternionCovariance inversePose3DQuaternionCovarianceQuaternion(
  const PoseQuaternionCovariance & pose_quaternion_covariance);

/*
  / @brief Inverse of a pose with covariance in RPY representation
 */
PoseRPYCovariance inversePose3DRPYCovarianceRPY(const PoseRPYCovariance & pose_rpy_covariance);

/*
  / @brief Inverse of a pose in quaternion with covariance in RPY representation (ROS convention)
 */
PoseQuaternionCovarianceRPY inversePose3DQuaternionCovarianceRPY(
  const PoseQuaternionCovarianceRPY & pose_quaternion_covariance_rpy);

}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP
