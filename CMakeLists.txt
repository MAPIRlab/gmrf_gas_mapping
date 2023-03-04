cmake_minimum_required(VERSION 3.5)
project(gmrf_gas_mapping)

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(olfaction_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_ros REQUIRED)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)


if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()


# let the compiler search for headers in the include folder
include_directories(include)

# Add Executables
add_executable(gmrf_node src/gmrf_map.cpp src/gmrf_node.cpp)
ament_target_dependencies(gmrf_node rclcpp std_msgs sensor_msgs geometry_msgs olfaction_msgs Eigen3 tf2 tf2_ros)

add_executable(gmrf_simple_publisher src/gmrf_simple_publisher.cpp)
ament_target_dependencies(gmrf_simple_publisher rclcpp std_msgs olfaction_msgs)


# Declare and install executables
install(TARGETS
  gmrf_node
  gmrf_simple_publisher
  DESTINATION lib/${PROJECT_NAME}
)

# Install launch files.
install(DIRECTORY
  launch
  maps
  DESTINATION share/${PROJECT_NAME}/
)

ament_package()