/*********************************************************************
*
* Software License Agreement (BSD License)
*
*  Copyright (c)  2015, Örebro University, Sweden
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.

* Author: Javier G. Monroy
* KernelDM+V Implementation:  Victor Hernandez
*********************************************************************/

//-----------------------------------------------------
// ROS2 wrapper for the GMRF gas-distribution mapping
//-----------------------------------------------------
#include <rclcpp/rclcpp.hpp>
#include "std_msgs/msg/float32.hpp"
#include "geometry_msgs/msg/pose_with_covariance_stamped.hpp"
#include "geometry_msgs/msg/transform_stamped.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "nav_msgs/msg/occupancy_grid.hpp"
#include "olfaction_msgs/msg/gas_sensor.hpp"
#include "olfaction_msgs/msg/observation.hpp"
#include "gmrf_gas_mapping/gmrf_map.h"

#include <boost/thread/mutex.hpp>
#include <boost/math/constants/constants.hpp>
#include <tf2_ros/transform_listener.h>
#include "tf2_ros/buffer.h"

/*
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/MarkerArray.h>
*/

//Services
#include "olfaction_msgs/srv/suggest_next_observation_location.h"



class Cgmrf: public rclcpp::Node
{
public:
    Cgmrf();
    ~Cgmrf();
    void publishMaps();

    //GMRF variables    
    CGMRF_map                               *my_map;			// The Online Gas Distribution Map being generated
    nav_msgs::msg::OccupancyGrid             occupancyMap;       //Occupancy GridMap of the environment

    //Node Params
    std::string                             sensor_topic, observation_topic;
    std::string                             frame_id;
    std::string                             occupancy_map_topic;
    double                                  cell_size;
    double                                  exec_freq;
    std::string                             colormap;
    int                                     max_pclpoints_cell;
    double                                  max_sensor_val;
    double                                  min_sensor_val;
    double                                  suggest_next_location_sensor_th;
    std::string                             output_csv_folder;        // Optional output CSV file with the latest map. Get continuously overwritten. Ignored if empty

    double                                  GMRF_lambdaPrior;       // [GMRF model] The information (Lambda) of prior factors
    double                                  GMRF_lambdaObs;         // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    double                                  GMRF_lambdaObsLoss;     // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)


    //Variables
    bool                                    module_init;
    boost::mutex                            mutex_nose;
    boost::mutex                            mutex_position;
    bool                                    new_data_position;
    float                                   curr_reading;
    float                                   curr_x;
    float                                   curr_y;



protected:
    //Subscriptions & Publishers
    rclcpp::Subscription<olfaction_msgs::msg::GasSensor>::SharedPtr sub_gas_sensor;
    rclcpp::Subscription<olfaction_msgs::msg::Observation>::SharedPtr sub_observation;
    rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr sub_occupancy_map;

    //Callbacks
    void sensorCallback(const olfaction_msgs::msg::GasSensor::SharedPtr msg);
    void mapCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg);
    void obsCallback(const olfaction_msgs::msg::Observation::SharedPtr msg);

    std::shared_ptr<tf2_ros::TransformListener> tf_listener_{nullptr};
    std::unique_ptr<tf2_ros::Buffer> tf_buffer_;
};



