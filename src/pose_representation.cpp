#include "covariance_geometry/pose_representation.hpp"

namespace covariance_geometry
{

  void Pose3DRPYTo3DQuaternion(const PoseRPY& pose_in, PoseQuaternion& pose_out)
  {
    // Copy position
    pose_out.first = pose_in.first;
    // Axis-angle to quaternion
    RPYToQuaternion(pose_in.second, pose_out.second);
  }

  void Pose3DQuaternionTo3DRPY(const PoseQuaternion& pose_in, PoseRPY& pose_out)
  {
    // Copy position
    pose_out.first = pose_in.first;
    // Quaternion to axis-angle
    QuaternionToRPY(pose_in.second, pose_out.second);
  }

  inline void RPYToQuaternion(const Eigen::Vector3d& rpy, Eigen::Quaterniond& quaternion)
  {
    // Axis-angle to quaternion
    Eigen::AngleAxisd rollAngle(rpy.x(), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(rpy.y(), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(rpy.z(), Eigen::Vector3d::UnitZ());
    quaternion = yawAngle * pitchAngle * rollAngle;
  }

  inline void QuaternionToRPY(const Eigen::Quaterniond& quaternion, Eigen::Vector3d& rpy)
  {
    // Quaternion to axis-angle
    rpy = quaternion.toRotationMatrix().eulerAngles(0, 1, 2);
  }

} // namespace covariance_geometry