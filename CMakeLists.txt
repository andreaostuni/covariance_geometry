cmake_minimum_required(VERSION 3.8)
project(covariance_geometry)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# set(CMAKE_BUILD_TYPE Debug)
find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)

set(dependencies
  Eigen3
  mrpt-poses
)

add_library(pose_covariance_geometry SHARED
  src/pose_covariance_composition.cpp
  src/pose_composition.cpp
  src/pose_representation.cpp
  src/pose_covariance_representation.cpp
  src/covariance_representation.cpp
)

include_directories(PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  # set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()

  find_package(ament_cmake_gtest REQUIRED)
  find_package(mrpt-poses REQUIRED)
  set(TEST_NAMES
    pose_composition_test
    pose_covariance_composition_test
    pose_representation_test
    pose_covariance_representation_test)
  foreach(TEST_NAME ${TEST_NAMES})
    ament_add_gtest(${TEST_NAME}
      test/${TEST_NAME}.cpp)
    ament_target_dependencies(${TEST_NAME}
      Eigen3
      mrpt-poses
    )
    target_link_libraries(${TEST_NAME} pose_covariance_geometry mrpt::poses)
  endforeach()
endif()

install(TARGETS pose_covariance_geometry
  EXPORT ${PROJECT_NAME}-export
  LIBRARY       DESTINATION lib
  PUBLIC_HEADER DESTINATION lib
)

install(DIRECTORY include/
  DESTINATION include
)

ament_export_include_directories("include/${PROJECT_NAME}" ${EIGEN3_INCLUDE_DIR})
ament_export_targets(${PROJECT_NAME}-export HAS_LIBRARY_TARGET)
ament_export_dependencies(${dependencies} Eigen3)

ament_package()
