cmake_minimum_required(VERSION 3.8)
project(covariance_geometry)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(CMAKE_BUILD_TYPE Debug)

find_package(Eigen3 REQUIRED)

add_library(pose_covariance_geometry SHARED 
  src/pose_covariance_composition.cpp
  src/pose_composition.cpp
  src/pose_representation.cpp
  src/pose_covariance_representation.cpp
  src/covariance_representation.cpp
  )

target_include_directories(pose_covariance_geometry PUBLIC include)

## tests
add_executable(test_node test/test.cpp)
target_include_directories(test_node PUBLIC include)
target_link_libraries(test_node pose_covariance_geometry)

# if(BUILD_TESTING)
#   find_package(ament_lint_auto REQUIRED)
#   # the following line skips the linter which checks for copyrights
#   # comment the line when a copyright and license is added to all source files
#   set(ament_cmake_copyright_FOUND TRUE)
#   # the following line skips cpplint (only works in a git repo)
#   # comment the line when this package is in a git repo and when
#   # a copyright and license is added to all source files
#   set(ament_cmake_cpplint_FOUND TRUE)
#   ament_lint_auto_find_test_dependencies()
# endif()

install(TARGETS pose_covariance_geometry EXPORT ${PROJECT_NAME}-export
  LIBRARY       DESTINATION lib
  PUBLIC_HEADER DESTINATION lib
)