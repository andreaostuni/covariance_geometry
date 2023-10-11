cmake_minimum_required(VERSION 3.15)
project(covariance_geometry_cpp)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(PROJECT_VERSION 0.0.1)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED)

include_directories(include)

add_library(${PROJECT_NAME} SHARED
  src/pose_representation.cpp
  src/covariance_representation.cpp
  src/pose_composition.cpp
  src/pose_covariance_representation.cpp
  src/pose_covariance_composition.cpp
)

target_include_directories(${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
)

if(BUILD_TESTING)
  find_package(GTest REQUIRED)
  find_package(mrpt-poses 2.10 REQUIRED)

  set(TEST_NAMES
  pose_composition_test
  pose_covariance_composition_test
  pose_representation_test
  pose_covariance_representation_test)

  include(GoogleTest)
  foreach(TEST_NAME ${TEST_NAMES})
    enable_testing()
    add_executable(${TEST_NAME}
      test/${TEST_NAME}.cpp)
    target_link_libraries(
      ${TEST_NAME}
      GTest::gtest_main
      ${PROJECT_NAME}
      mrpt::poses
    )
    include(GoogleTest)
    gtest_discover_tests(${TEST_NAME})
  endforeach()
endif()

include(GNUInstallDirs)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION include)

export(TARGETS ${PROJECT_NAME}
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake)
set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE ${PROJECT_NAME})

install(EXPORT ${PROJECT_NAME}-targets
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

include(CMakePackageConfigHelpers)
# generate the config file that includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
