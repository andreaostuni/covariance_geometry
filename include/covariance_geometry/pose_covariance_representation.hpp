#ifndef COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP
#define COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP

#include "covariance_geometry/pose_representation.hpp"
#include "covariance_geometry/covariance_representation.hpp"

namespace covariance_geometry{
  using PoseQuaternionCovariance = std::pair<PoseQuaternion, Eigen::Matrix7d>;
  using PoseRPYCovariance = std::pair<PoseRPY, Eigen::Matrix6d>;
  using PoseQuaternionCovarianceRPY = std::pair<PoseQuaternion, Eigen::Matrix6d>; // Ros convention

  /*
  / @brief Convert a Pose with covariance matrix from quaternion representation to RPY representation
  */
  void Pose3DQuaternionCovarianceTo3DRPYCovariance(const PoseQuaternionCovariance& pose_quaternion_covariance, PoseRPYCovariance& pose_rpy_covariance);

  /*
  / @brief Convert a Pose with covariance matrix from RPY representation to quaternion representation
  */
  void Pose3DRPYCovarianceTo3DQuaternionCovariance(const PoseRPYCovariance& pose_rpy_covariance, PoseQuaternionCovariance& pose_quaternion_covariance);

  /*
  / @brief Convert a Pose with covariance matrix from RPY representation to quaternion representation but covariance still in RPY (ROS convention)
  */
  void Pose3DRPYCovarianceTo3DQuaternionCovarianceRPY(const PoseRPYCovariance& pose_rpy_covariance, PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy);

  /*
  / @brief Convert a Pose with covariance matrix from quaternion representation to quaternion representation but covariance still in RPY (ROS convention)
  */
  void Pose3DQuaternionCovarianceTo3DQuaternionCovarianceRPY(const PoseQuaternionCovariance& pose_quaternion_covariance, PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy);

  /*
  / @brief Convert a Pose with covariance matrix from quaternion representation with covariance in RPY (ROS representation) to RPY representation
  */
  void Pose3DQuaternionCovarianceRPYTo3DRPYCovariance(const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy, PoseRPYCovariance& pose_rpy_covariance);

  /*
  / @brief Convert a Pose with covariance matrix from quaternion representation with covariance in RPY (ROS representation) to quaternion representation 
  */
  void Pose3DQuaternionCovarianceRPYTo3DQuaternionCovariance(const PoseQuaternionCovarianceRPY& pose_quaternion_covariance_rpy, PoseQuaternionCovariance& pose_quaternion_covariance);
  
} // namespace covariance_geometry

#endif // COVARIANCE_GEOMETRY_POSE_COVARIANCE_REPRESENTATION_HPP