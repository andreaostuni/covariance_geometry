#include "covariance_geometry/pose_covariance_representation.hpp"

namespace covariance_geometry{

  void Pose3DQuaternionCovarianceTo3DRPYCovariance(const PoseQuaternionCovariance& pose_quaternion_covariance, PoseRPYCovariance& pose_rpy_covariance){
    // Convert pose
    Pose3DQuaternionTo3DRPY(pose_quaternion_covariance.first, pose_rpy_covariance.first);

    // Convert covariance
    covariance3DQuaternionTo3DRPY(pose_quaternion_covariance.first.second, pose_quaternion_covariance.second, pose_rpy_covariance.second);
  }

  void Pose3DRPYCovarianceTo3DQuaternionCovariance(const PoseRPYCovariance& pose_rpy_covariance, PoseQuaternionCovariance& pose_quaternion_covariance){
    // Convert pose
    Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance.first);

    // Convert covariance
    covariance3DRPYTo3DQuaternion(pose_rpy_covariance.first.second, pose_rpy_covariance.second, pose_quaternion_covariance.second);
  }

  void Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(const PoseRPYCovariance& pose_rpy_covariance, PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy){
    // Convert pose
    Pose3DRPYTo3DQuaternion(pose_rpy_covariance.first, pose_quaternion_covariance_rpy.first);
  }

  void Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(const PoseQuaternionCovariance& pose_quaternion_covariance, PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy){
    // Convert covariance
    covariance3DQuaternionTo3DRPY(pose_quaternion_covariance.first.second, pose_quaternion_covariance.second, pose_quaternion_covariance_rpy.second);
  }

  void Pose3DQuaternionCovarianceRPYTo3DRPYCovariance(const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy, PoseRPYCovariance& pose_rpy_covariance){
    // Convert pose
    Pose3DQuaternionTo3DRPY(pose_quaternion_covariance_rpy.first, pose_rpy_covariance.first);
  }

  void Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy, PoseQuaternionCovariance& pose_quaternion_covariance){
    // Convert pose from quat to RPY
    Eigen::Vector3d pose_rpy;
    QuaternionToRPY(pose_quaternion_covariance_rpy.first.second, pose_rpy);
    // Convert covariance
    covariance3DRPYTo3DQuaternion(pose_rpy, pose_quaternion_covariance_rpy.second, pose_quaternion_covariance.second);
  }

} // namespace covariance_geometry

