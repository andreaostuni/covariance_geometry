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

void ComposePose3DQuaternion(
  const PoseQuaternion & a, const PoseQuaternion & b, PoseQuaternion & pose_out);

/*
  / @brief Pose composition for 6D poses in RPY form
  */

void ComposePose3DRPY(const PoseRPY & a, const PoseRPY & b, PoseRPY & pose_out);

}  // namespace covariance_geometry

#endif  // COVARIANCE_GEOMETRY_POSE_COMPOSITION_HPP
