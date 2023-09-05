#include <mrpt/poses/CPose3D.h>
#include <mrpt/poses/CPose3DPDFGaussian.h>
#include <mrpt/poses/CPose3DQuat.h>

#include <cmath>
#include <iostream>

#include "covariance_geometry/pose_composition.hpp"
#include "covariance_geometry/pose_covariance_composition.hpp"
#include "covariance_geometry/pose_covariance_representation.hpp"
#include "covariance_geometry/pose_representation.hpp"
#include "gtest/gtest.h"

using namespace covariance_geometry;
#define GTEST_COUT std::cerr << "[          ] [ INFO ]"

template<typename T>
bool isApprox(const T & a, const T & b, double epsilon = 1e-6)
{
  if constexpr (std::is_same_v<T, PoseQuaternion>) {
    GTEST_COUT << "a.first: " << a.first.transpose() << "\t";
    std::cout << "a.second: " << a.second << std::endl;
    GTEST_COUT << "b.first: " << b.first.transpose() << "\t";
    std::cout << "b.second: " << b.second << std::endl;
  } else {
    GTEST_COUT << "a.first: " << a.first.transpose() << "\t";
    std::cout << "a.second: " << a.second.transpose() << std::endl;
    GTEST_COUT << "b.first: " << b.first.transpose() << "\t";
    std::cout << "b.second: " << b.second.transpose() << std::endl;
  }
  return a.first.isApprox(b.first, epsilon) && a.second.isApprox(b.second, epsilon);
}

template<typename T>
bool MRPTtoEigen(const mrpt::poses::CPose3D & mrpt_pose, T & eigen_pose)
{
  if constexpr (std::is_same_v<T, PoseRPY>) {
    // eigen_pose.first = mrpt_pose.m_coords.asEigen();
    // mrpt_pose.getYawPitchRoll(eigen_pose.second(2), eigen_pose.second(1), eigen_pose.second(0));
    eigen_pose.first.x() = mrpt_pose.x();
    eigen_pose.first.y() = mrpt_pose.y();
    eigen_pose.first.z() = mrpt_pose.z();
    eigen_pose.second.x() = mrpt_pose.roll();
    eigen_pose.second.y() = mrpt_pose.pitch();
    eigen_pose.second.z() = mrpt_pose.yaw();
    return true;
  } else if constexpr (std::is_same_v<T, PoseQuaternion>) {
    eigen_pose.first = mrpt_pose.m_coords.asEigen();
    mrpt::math::CQuaternionDouble q;
    mrpt_pose.getAsQuaternion(q);
    eigen_pose.second = Eigen::Quaterniond(q.r(), q.x(), q.y(), q.z());
    return true;
  }
  throw std::runtime_error("MRPTtoEigen: Not compatible types");
  return false;
}

template<typename T>
bool EigenToMRPT(const T & eigen_pose, mrpt::poses::CPose3D & mrpt_pose)
{
  if constexpr (std::is_same_v<T, PoseRPY>) {
    mrpt_pose = mrpt_pose.FromXYZYawPitchRoll(
      eigen_pose.first(0), eigen_pose.first(1), eigen_pose.first(2), eigen_pose.second(2),
      eigen_pose.second(1), eigen_pose.second(0));
    return true;
  } else if constexpr (std::is_same_v<T, PoseQuaternion>) {
    mrpt::math::CQuaternionDouble q{
      eigen_pose.second.w(), eigen_pose.second.x(), eigen_pose.second.y(), eigen_pose.second.z()};
    mrpt_pose = mrpt_pose.FromQuaternionAndTranslation(
      q, eigen_pose.first.x(), eigen_pose.first.y(), eigen_pose.first.z());
    return true;
  }

  throw std::runtime_error("EigenToMRPT: Not compatible types");
  return false;
}

TEST(EigenMRPT, EigenToMRPTRPY)
{
  PoseRPY pr1, pr_out;
  mrpt::poses::CPose3D mrpt_pose;
  pr1.first = {1, 0, 0};
  pr1.second = {0.9128709, 0.4082483, 0};
  EigenToMRPT(pr1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pr_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, pr_out));

  pr1.first = {1, 0, 0};
  pr1.second = {0.9067432, 0.4055079, 0.1055943};
  EigenToMRPT(pr1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pr_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, pr_out));
}

TEST(EigenMRPT, MRPTtoEigenQuat)
{
  PoseQuaternion pq1, pq_out;
  mrpt::poses::CPose3D mrpt_pose;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  EigenToMRPT(pq1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pq_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  EigenToMRPT(pq1, mrpt_pose);
  MRPTtoEigen(mrpt_pose, pq_out);
  GTEST_COUT << "mrpt: " << mrpt_pose.asString() << std::endl;
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));
}

TEST(PoseComposition, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1, identity, p_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  identity.first = {0, 0, 0};
  identity.second = {1, 0, 0, 0};  // W X Y Z
  ComposePose3DQuaternion(pq1, identity, p_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, p_out));
}

TEST(PoseComposition, HandleZeroPoseRPY)
{
  PoseRPY pr1, identity, p_out;
  pr1.first = {1.0, 0.0, 0.0};
  pr1.second = {0.9128709, 0.4082483, 0};
  identity.first = {0.0, 0.0, 0.0};
  identity.second = {0.0, 0.0, 0.0};
  ComposePose3DRPY(pr1, identity, p_out);
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, p_out));
}

TEST(PoseComposition, HandleCompositionPoseQuaternion)
{
  PoseQuaternion pq1, pq2, pq3, p_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0.0, 0.0};
  pq2.first = {2, 0, 0};
  pq2.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  ComposePose3DQuaternion(pq1, pq2, p_out);

  mrpt::poses::CPose3D p1, p2;
  EigenToMRPT(pq1, p1);
  EigenToMRPT(pq2, p2);
  mrpt::poses::CPose3D p3 = p1 + p2;
  MRPTtoEigen(p3, pq3);

  EXPECT_TRUE(isApprox<PoseQuaternion>(pq3, p_out));
}

TEST(PoseComposition, HandleCompositionPoseRPY)
{
  PoseRPY pr1, pr2, pr3, p_out;
  pr1.first = {1, 0, 0};
  // pr1.second = {0.9128709, 0.4082483, 0};
  pr1.second = {1.5709, 0.0, 0.0};

  pr2.first = {2, 0, 0};
  // pr2.second = {0.9128709, 0.4082483, 0};
  pr2.second = {0.0, 1.5709, 0.0};

  ComposePose3DRPY(pr1, pr2, p_out);

  mrpt::poses::CPose3D p1, p2, p3;
  // mrpt::poses::CPose3D p4, p5, p6;
  // pr1.first = {1, 0, 0};
  // pr1.second = {0.9128709, 0.4082483, 0};
  // pr2.first = {2, 0, 0};
  // pr2.second = {0.9128709, 0.4082483, 0};

  EigenToMRPT(pr1, p1);
  EigenToMRPT(pr2, p2);
  p3 = p1 + p2;
  MRPTtoEigen(p3, pr3);

  EXPECT_TRUE(isApprox<PoseRPY>(pr3, p_out));

  // pr1.first = {1, 0, 0};
  // pr1.second = {0.9128709, 0.4082483, 0};
  // pr2.first = {2, 0, 0};
  // pr2.second = {0.9067432, 0.4055079, 0.1055943};

  // ComposePose3DRPY(pr1, pr2, p_out);

  // EigenToMRPT(pr1, p4);
  // EigenToMRPT(pr2, p5);
  // p6 = p4 + p5;
  // MRPTtoEigen(p6, pr3);

  // EXPECT_TRUE(isApprox<PoseRPY>(pr3, p_out));
}

TEST(PoseConversion, HandleZeroPoseQuaternion)
{
  PoseQuaternion pq1;
  PoseRPY pr1, p_out;
  pq1.first = {0, 0, 0};
  pq1.second = {1, 0, 0, 0};
  pr1.first = {0, 0, 0};
  pr1.second = {0, 0, 0};
  Pose3DQuaternionTo3DRPY(pq1, p_out);
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, p_out));
}

TEST(PoseConversion, HandleZeroPoseRPY)
{
  PoseQuaternion pq1, p_out;
  PoseRPY pr1;
  pq1.first = {0, 0, 0};
  pq1.second = {1, 0, 0, 0};
  pr1.first = {0, 0, 0};
  pr1.second = {0, 0, 0};
  Pose3DRPYTo3DQuaternion(pr1, p_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, p_out));
}

TEST(PoseConversion, HandleConversionPoseQuaternion)
{
  PoseQuaternion pq1;
  PoseRPY pr1, pr_out;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8410687, 0.0, 0};
  Pose3DQuaternionTo3DRPY(pq1, pr_out);
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, pr_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8545253, 0.1538007, 0.1742029};
  Pose3DQuaternionTo3DRPY(pq1, pr_out);
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, pr_out));
}

TEST(PoseConversion, HandleConversionPoseRPY)
{
  PoseQuaternion pq1, pq_out;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8410687, 0.0, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.9067432, 0.4055079, 0.1055943, 0.0472232};
  pr1.first = {1, 0, 0};
  pr1.second = {0.8545253, 0.1538007, 0.1742029};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));
}

TEST(PoseConversion, HandleConversionPoseRPY_90Pitch)
{
  PoseQuaternion pq1, pq_out;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.7071068, 0.0, 0.7071068, 0.0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.0, M_PI_2, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));

  pq1.first = {1, 0, 0};
  pq1.second = {0.7071068, 0.0, -0.7071068, 0.0};
  pr1.first = {1, 0, 0};
  pr1.second = {0.0, -M_PI_2, 0.0};
  Pose3DRPYTo3DQuaternion(pr1, pq_out);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq_out));
}

TEST(PoseConversion, CyclicPoseConversionQuaternion)
{
  PoseQuaternion pq1, pq2;
  PoseRPY pr1;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  Pose3DQuaternionTo3DRPY(pq1, pr1);
  Pose3DRPYTo3DQuaternion(pr1, pq2);
  EXPECT_TRUE(isApprox<PoseQuaternion>(pq1, pq2));
}

TEST(PoseConversion, CyclicPoseConversionRPY)
{
  PoseQuaternion pq1;
  PoseRPY pr1, pr2;
  pq1.first = {1, 0, 0};
  pq1.second = {0.9128709, 0.4082483, 0, 0};
  Pose3DRPYTo3DQuaternion(pr1, pq1);
  Pose3DQuaternionTo3DRPY(pq1, pr2);
  EXPECT_TRUE(isApprox<PoseRPY>(pr1, pr2));
}

void test_pose_composition_conversion()
{
  std::cout << "Testing pose composition and conversion.." << std::endl;
  PoseQuaternion pq1, pq2, pq3;
  PoseRPY pr1, pr2, pr3;
  Eigen::Vector4d q1(0, 0.4082483, 0, 0.9128709);
  //!
  //! if initialized as Eigen::Quaterniond q1(0.9128709, 0.4082483, 0, 0);, the
  //! order is w x y z
  //! if initialized as Eigen::Quaterniond = Eigen::Vector4d q1(0.9128709,
  //! 0.4082483, 0, 0);, the order is x y z w
  pq1.first = {1, 0, 0};
  pq1.second = q1;

  pq2.first = {2, 0, 0};
  pq2.second = q1;
  ComposePose3DQuaternion(pq1, pq2, pq3);

  std::cout << "Pose1 t: "
            << "\n"
            << pq1.first << "\n"
            << "q: " << pq1.second << std::endl;
  std::cout << "Pose2 t: "
            << "\n"
            << pq2.first << "\n"
            << "q: " << pq2.second << std::endl;
  std::cout << "Composition t: "
            << "\n"
            << pq3.first << "\n"
            << "q: " << pq3.second << std::endl;

  Pose3DQuaternionTo3DRPY(pq1, pr1);
  Pose3DQuaternionTo3DRPY(pq2, pr2);
  ComposePose3DRPY(pr1, pr2, pr3);

  std::cout << "Pose1 t: "
            << "\n"
            << pr1.first << "\n"
            << "rpy: " << pr1.second << std::endl;
  std::cout << "Pose2 t: "
            << "\n"
            << pr2.first << "\n"
            << "rpy: " << pr2.second << std::endl;
  std::cout << "Composition t: "
            << "\n"
            << pr3.first << "\n"
            << "rpy: " << pr3.second << std::endl;

  std::cout << "Converting back.." << std::endl;
  PoseQuaternion pq4;
  Pose3DRPYTo3DQuaternion(pr3, pq4);
  std::cout << "Composition t: "
            << "\n"
            << pq3.first << "\n"
            << "q: " << pq3.second << std::endl;
  std::cout << "Composition t: "
            << "\n"
            << pq4.first << "\n"
            << "q: " << pq4.second << std::endl;
}

void test_pose_with_covariance_composition_conversion()
{
  std::cout << "Testing pose with covariance composition and conversion.." << std::endl;

  // ROS type
  PoseQuaternionCovarianceRPY p1, p2, p3;
  // Eigen::Vector4d q1(0, 0.4082483, 0, 0.9128709);
  Eigen::Vector4d q1(0.0, 0.0, 0.0, 1.0);
  Eigen::Vector3d p1_t(1.0, 0.0, 0.0);
  Eigen::Quaterniond p1_q;
  p1_q = q1;
  auto p1_pose = std::make_pair(p1_t, p1_q);
  p1 = std::make_pair(p1_pose, Eigen::DiagonalMatrix<double, 6, 6>(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));

  Eigen::Vector3d p2_t(1.0, 0.0, 0.0);
  Eigen::Quaterniond p2_q;
  p2_q = q1;
  auto p2_pose = std::make_pair(p2_t, p2_q);
  p2 = std::make_pair(p2_pose, Eigen::DiagonalMatrix<double, 6, 6>(1.0, 1.0, 1.0, 1.0, 1.0, 1.0));

  ComposePoseQuaternionCovarianceRPY(p1, p2, p3);

  std::cout << "Pose1 t: "
            << "\n"
            << p1.first.first << "\n"
            << "q: \n"
            << p1.first.second
            << "\n"
    "cov1: \n"
            << p1.second << std::endl;

  std::cout << "Pose2 t: "
            << "\n"
            << p2.first.first << "\n"
            << "q: \n"
            << p2.first.second
            << "\n"
    "cov2: \n"
            << p2.second << std::endl;

  std::cout << "Composition t: "
            << "\n"
            << p3.first.first << "\n"
            << "q: \n"
            << p3.first.second
            << "\n"
    "cov: \n"
            << p3.second << std::endl;
}

int main(int argc, char ** argv)
{
  // test_pose_composition_conversion();
  // test_pose_with_covariance_composition_conversion();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  // return 0;
}
