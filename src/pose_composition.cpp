#include "covariance_geometry/pose_composition.hpp"
#include "covariance_geometry/pose_representation.hpp"

namespace covariance_geometry
{

  void ComposePose3DQuaternion(const PoseQuaternion &a, const PoseQuaternion &b, PoseQuaternion &pose_out)
  {
    // Position composition
    pose_out.first = a.first + a.second * b.first;
    // Quaternion composition
    pose_out.second = a.second * b.second;
  }

  void ComposePose3DRPY(const PoseRPY &a, const PoseRPY &b, PoseRPY &pose_out)
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

} // namespace covariance_geometry