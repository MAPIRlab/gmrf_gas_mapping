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
#cv_bridge
#pcl_ros
#tf  
#visualization_msgs  


# let the compiler search for headers in the include folder
include_directories(include)

# Add Executables
add_executable(gmrf_gas_mapping_node src/gmrf_map.cpp src/gmrf_node.cpp)
ament_target_dependencies(gmrf_gas_mapping_node rclcpp std_msgs sensor_msgs geometry_msgs olfaction_msgs Eigen3 tf2 tf2_ros)


# Declare and install executables
install(TARGETS
  gmrf_gas_mapping_node
  DESTINATION lib/${PROJECT_NAME}
)

# Install launch files.
install(DIRECTORY
  launch
  DESTINATION share/${PROJECT_NAME}/
)

ament_package()
