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

void Pose3DQuaternionCovarianceTo3DRPYCovariance(const PoseQuaternionCovariance& pose_quaternion_covariance,
                                                 PoseRPYCovariance& pose_rpy_covariance)
{
  // Convert pose
  Pose3DQuaternionTo3DRPY(pose_quaternion_covariance.first, pose_rpy_covariance.first);

  // Convert covariance
  covariance3DQuaternionTo3DRPY(pose_quaternion_covariance.first.second, pose_quaternion_covariance.second,
                                pose_rpy_covariance.second);
}

void Pose3DRPYCovarianceTo3DQuaternionCovariance(const PoseRPYCovariance& pose_rpy_covariance,
                                                 PoseQuaternionCovariance& pose_quaternion_covariance)
{
  // Convert pose
  Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance.first);

  // Convert covariance
  covariance3DRPYTo3DQuaternion(pose_rpy_covariance.first.second, pose_rpy_covariance.second,
                                pose_quaternion_covariance.second);
}

void Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(const PoseRPYCovariance& pose_rpy_covariance,
                                                    PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy)
{
  // Convert pose
  Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance_rpy.first);
}

void Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(const PoseQuaternionCovariance& pose_quaternion_covariance,
                                                           PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy)
{
  // Copy pose
  pose_quaternion_covariance_rpy.first = pose_quaternion_covariance.first;
  // Convert covariance
  covariance3DQuaternionTo3DRPY(pose_quaternion_covariance.first.second, pose_quaternion_covariance.second,
                                pose_quaternion_covariance_rpy.second);
}

void Pose3DQuaternionCovarianceRPYTo3DRPYCovariance(const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy,
                                                    PoseRPYCovariance& pose_rpy_covariance)
{
  // Convert pose
  Pose3DQuaternionTo3DRPY(pose_quaternion_covariance_rpy.first, pose_rpy_covariance.first);
}

void Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(
    const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy,
    PoseQuaternionCovariance& pose_quaternion_covariance)
{
  // TODO: check this method
  // Copy pose
  pose_quaternion_covariance.first = pose_quaternion_covariance_rpy.first;

  // Convert pose rotation from quat to RPY
  Eigen::Vector3d pose_rpy;
  QuaternionToRPY(pose_quaternion_covariance_rpy.first.second, pose_rpy);
  // Convert covariance
  covariance3DRPYTo3DQuaternion(pose_rpy, pose_quaternion_covariance_rpy.second, pose_quaternion_covariance.second);
}

}  // namespace covariance_geometry
