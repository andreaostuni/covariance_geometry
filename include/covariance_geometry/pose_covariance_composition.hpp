#ifndef COVARIANCE_GEOMETRY_POSE_COVARIANCE_COMPOSITION_HPP
#define COVARIANCE_GEOMETRY_POSE_COVARIANCE_COMPOSITION_HPP

#include <utility>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace Eigen{
  using Matrix6d = Eigen::Matrix<double, 6, 6>;
  using Matrix7d = Eigen::Matrix<double, 7, 7>;
  using Matrix3_7d = Eigen::Matrix<double, 3, 7>;
  using Matrix3_4d = Eigen::Matrix<double, 3, 4>;
}

namespace covariance_geometry
{

    using PoseQuaternion = std::pair<Eigen::Vector3d, Eigen::Quaterniond>;
    using PoseRPY = std::pair<Eigen::Vector3d, Eigen::Vector3d>;
    using PoseQuaternionCovariance = std::pair<PoseQuaternion, Eigen::Matrix7d>;
    using PoseRPYCovariance = std::pair<PoseRPY, Eigen::Matrix6d>;
    using PoseQuaternionCovarianceRPY = std::pair<PoseQuaternion, Eigen::Matrix6d>; // Ros convention

  /*
  / @brief Pose with covariance composition for 6D poses with 6x6 covariance
  */
  void ComposePoseRPYCovarianceRPY(const PoseRPYCovariance &a, const PoseRPYCovariance &b, PoseRPYCovariance &out);
  
  /*
  / @brief Pose with covariance composition for 7D poses with 6x6 covariance (ROS convention)
  */
  void ComposePoseQuaternionCovarianceRPY(const PoseQuaternionCovarianceRPY &a, const PoseQuaternionCovarianceRPY &b, PoseQuaternionCovarianceRPY &out);

  /*
  / @brief Pose with covariance composition for 7D poses with 7x7 covariance
  */
  void ComposePoseQuaternionCovariance(const PoseQuaternionCovariance &a, const PoseQuaternionCovariance &b, PoseQuaternionCovariance &out);

  /*
  / @brief Compute pose-pose composition function jacobian wrt pose_a
  */
  void JacobianPosePoseCompositionA(const PoseQuaternion &pose_a, const PoseQuaternion &pose_b, Eigen::Matrix7d &jacobian);

  /*
  / @brief Compute pose-pose composition function jacobian wrt pose_b
  */
  void JacobianPosePoseCompositionB(const PoseQuaternion &pose_a, Eigen::Matrix7d &jacobian);

  /*
  / @brief Compute pose-point composition function jacobian wrt pose
  */
  void JacobianPosePointComposition(const PoseQuaternion &pose, const Eigen::Vector3d &point, Eigen::Matrix3_7d &jacobian);

  /*
  / @brief Compute pose-point composition function jacobian wrt point
  */
  void JacobianPosePointComposition(const PoseQuaternion &pose, Eigen::Matrix3d &jacobian);

  /*
  / @brief Compute pose-point composition function jacobian wrt pose quaternion
  */
  void JacobianQuaternionPointComposition(const Eigen::Quaterniond &quaternion, const Eigen::Vector3d &point, Eigen::Matrix3_4d &jacobian);

} // namespace covariance_geometry

#endif // COVARIANCE_GEOMETRY_POSE_COVARIANCE_COMPOSITION_HPP