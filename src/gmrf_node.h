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
// ROS wrapper for the GMRF gas-distribution mapping
//-----------------------------------------------------
#include "ros/ros.h"
#include "std_msgs/Float32.h"
#include "nav_msgs/Odometry.h"
#include "nav_msgs/OccupancyGrid.h"
#include <sensor_msgs/PointCloud2.h>
#include "olfaction_msgs/gas_sensor.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include <visualization_msgs/MarkerArray.h>
#include <boost/thread/mutex.hpp>
#include <boost/math/constants/constants.hpp>
#include <tf/transform_listener.h>
#include "std_msgs/Float32MultiArray.h"

#include "gmrf_map.h"

//Services
#include "olfaction_msgs/suggestNextObservationLocation.h"



class Cgmrf
{
public:
    Cgmrf();
    ~Cgmrf();
    void publishMaps();

    //GMRF variables    
    CGMRF_map                               *my_map;			// The Online Gas Distribution Map being generated
    nav_msgs::OccupancyGrid                 occupancyMap;       //Occupancy GridMap of the environment

    //Node Params
    std::string                             sensor_topic;
    std::string                             frame_id;
    std::string                             occupancy_map_topic;
    double                                  cell_size;
    double                                  exec_freq;
    std::string                             colormap;
    int                                     max_pclpoints_cell;
    double                                  max_sensor_val;
    double                                  min_sensor_val;
    double                                  suggest_next_location_sensor_th;
    std::string                             output_csv_file;        // Optional output CSV file with the latest map. Get continuously overwritten. Ignored if empty

    double                                  GMRF_lambdaPrior;       // [GMRF model] The information (Lambda) of prior factors
    double                                  GMRF_lambdaObs;         // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    double                                  GMRF_lambdaObsLoss;     // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)


    //Variables
    bool module_init;
    boost::mutex                            mutex_nose;
    boost::mutex                            mutex_position;
    bool                                    new_data_position;
    float                                   curr_reading;
    float                                   curr_x;
    float                                   curr_y;



protected:
    ros::NodeHandle                         n;
    tf::TransformListener                   tf_listener; //Do not put inside the callback

    //Subscriptions & Publishers
    ros::Subscriber sub_sensor, occupancyMap_sub;
    ros::Publisher mean_advertise, var_advertise, mean_raw_advertise, var_raw_advertise;

    //Callbacks
    void sensorCallback(const olfaction_msgs::gas_sensorPtr msg);
    void mapCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg);
};



