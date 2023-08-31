#include "covariance_geometry/pose_covariance_composition.hpp"
#include "covariance_geometry/pose_composition.hpp"
#include "covariance_geometry/pose_representation.hpp"

#include <iostream>

using namespace covariance_geometry;


void test_pose_composition_conversion()
{
    std::cout << "Testing pose composition and conversion.." << std::endl;
    PoseQuaternion pq1, pq2, pq3;
    PoseRPY pr1, pr2, pr3;

    pq1.first = {1, 0, 0};
    pq1.second = {0.9128709, 0.4082483, 0, 0};
    pq2.first = {2, 0, 0};
    pq2.second = {0.9128709, 0.4082483, 0, 0};
    ComposePose3DQuaternion(pq1, pq2, pq3);

    std::cout << "Pose1 t: "       << "\n" << pq1.first << "\n" << "q: " << pq1.second << std::endl;
    std::cout << "Pose2 t: "       << "\n" << pq2.first << "\n" << "q: " << pq2.second << std::endl;
    std::cout << "Composition t: " << "\n" << pq3.first << "\n" << "q: " << pq3.second << std::endl;

    Pose3DQuaternionTo3DRPY(pq1, pr1);
    Pose3DQuaternionTo3DRPY(pq2, pr2);
    ComposePose3DRPY(pr1, pr2, pr3);

    std::cout << "Pose1 t: "       << "\n" << pr1.first << "\n" << "rpy: " << pr1.second << std::endl;
    std::cout << "Pose2 t: "       << "\n" << pr2.first << "\n" << "rpy: " << pr2.second << std::endl;
    std::cout << "Composition t: " << "\n" << pr3.first << "\n" << "rpy: " << pr3.second << std::endl;


    std::cout << "Converting back.." << std::endl;
    PoseQuaternion pq4;
    Pose3DRPYTo3DQuaternion(pr3, pq4);
    std::cout << "Composition t: " << "\n" << pq3.first << "\n" << "q: " << pq3.second << std::endl;
    std::cout << "Composition t: " << "\n" << pq4.first << "\n" << "q: " << pq4.second << std::endl;
}

void test_pose_with_covariance_composition_conversion()
{
    std::cout << "Testing pose with covariance composition and conversion.." << std::endl;

    // ROS type
    PoseQuaternionCovarianceRPY p1, p2, p3;
    p1.first.first = {1, 0, 0};
    p1.first.second = {0.9128709, 0.4082483, 0, 0};
    p1.second = Eigen::DiagonalMatrix<double, 6, 6>::DiagonalMatrix(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

    p2.first.first = {2, 0, 0};
    p2.first.second = {0.9128709, 0.4082483, 0, 0};
    p2.second = Eigen::DiagonalMatrix<double, 6, 6>::DiagonalMatrix(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);    

    ComposePoseQuaternionCovarianceRPY(p1, p2, p3);

    std::cout << "Pose1 t: "       << "\n" << /
    p1.first.first << "\n" << "q: " << p1.first.second << /
    pstd::endl;
    std::cout << "Pose2 t: "       << "\n" << pq2.first << "\n" << "q: " << pq2.second << std::endl;
    std::cout << "Composition t: " << "\n" << pq3.first << "\n" << "q: " << pq3.second << std::endl;




}

int main()
{
    test_pose_composition_conversion();
    return 0;
}