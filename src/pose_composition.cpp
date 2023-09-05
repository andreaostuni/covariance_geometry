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

#include "covariance_geometry/pose_composition.hpp"

#include "covariance_geometry/pose_representation.hpp"

namespace covariance_geometry
{

void ComposePose3DQuaternion(const PoseQuaternion& a, const PoseQuaternion& b, PoseQuaternion& pose_out)
{
  // Position composition
  pose_out.first = a.first + a.second * b.first;
  // Quaternion composition
  pose_out.second = a.second * b.second;
}

void ComposePose3DRPY(const PoseRPY& a, const PoseRPY& b, PoseRPY& pose_out)
{
  // Convert to quaternion
  PoseQuaternion a_quaternion, b_quaternion, pose_out_quaternion;
  Pose3DRPYTo3DQuaternion(a, a_quaternion);
  Pose3DRPYTo3DQuaternion(b, b_quaternion);
  // Compose
  ComposePose3DQuaternion(a_quaternion, b_quaternion, pose_out_quaternion);
  // Convert back to RPY
  Pose3DQuaternionTo3DRPY(pose_out_quaternion, pose_out);
}

}  // namespace covariance_geometry
