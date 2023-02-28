//========================================================================================
//	GMRF_node
//	Description: implements the Gaussian Markov rando field gas-distribution mapping algorithm.
//  See: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/193	
//		
//	topics subscribed:
//  topics published:
//	services:
//				
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//	Revision log:
//	version: 1.0	16/09/2016
//========================================================================================

#include "gmrf_gas_mapping/gmrf_node.h"

using std::placeholders::_1;

// -------------
// Cgmrf
//--------------
Cgmrf::Cgmrf(): Node("GMRF_node")
{    
    RCLCPP_INFO(this->get_logger(), "[GMRF] Loading parameters");

    //-----------------------------
    // Declare and Load Parameters
    //-----------------------------    
    frame_id = this->declare_parameter<std::string>("frame_id", "map");
    occupancy_map_topic = this->declare_parameter<std::string>("occupancy_map_topic", "/map");
    sensor_topic = this->declare_parameter<std::string>("sensor_topic", "/Mox01/Sensor_reading");
    output_csv_folder = this->declare_parameter<std::string>("output_csv_folder", "");
    exec_freq = this->declare_parameter<double>("exec_freq", 2.0);
    cell_size = this->declare_parameter<double>("cell_size", 0.5);
    min_sensor_val = this->declare_parameter<double>("min_sensor_val", 0.0);
    max_sensor_val = this->declare_parameter<double>("max_sensor_val", 5.0);

    GMRF_lambdaPrior = this->declare_parameter<double>("GMRF_lambdaPrior", 0.5); // [GMRF model] The information (Lambda) of prior factors
    GMRF_lambdaObs = this->declare_parameter<double>("GMRF_lambdaObs", 10.0);   // [GMRF model] The initial information (Lambda) of each observation (this information will decrease with time)
    GMRF_lambdaObsLoss = this->declare_parameter<double>("GMRF_lambdaObsLoss", 0.0);    // [GMRF model] The loss of information (Lambda) of the observations with each iteration (see AppTick)
    //param_n.param<double>("suggest_next_location_sensor_th", suggest_next_location_sensor_th, 0.1);

    // Display main parameters
    RCLCPP_INFO(this->get_logger(), "[GMRF] frame_id: %s", frame_id.c_str());
    RCLCPP_INFO(this->get_logger(), "[GMRF] occupancy_map_topic: %s",  occupancy_map_topic.c_str());
    RCLCPP_INFO(this->get_logger(), "[GMRF] sensor_topic: %s",  sensor_topic.c_str());
    RCLCPP_INFO(this->get_logger(), "[GMRF] output_csv_folder: %s",  output_csv_folder.c_str());

    //----------------------------------
    // Subscriptions
    //----------------------------------
    sub_gas_sensor = this->create_subscription<olfaction_msgs::msg::GasSensor>(sensor_topic ,1, std::bind(&Cgmrf::sensorCallback, this, _1) );
    sub_occupancy_map = this->create_subscription<nav_msgs::msg::OccupancyGrid>(occupancy_map_topic,1, std::bind(&Cgmrf::mapCallback, this, _1) );
 

    tf_buffer_ = std::make_unique<tf2_ros::Buffer>(this->get_clock());
    tf_listener_ = std::make_shared<tf2_ros::TransformListener>(*tf_buffer_);

    //----------------------------------
    // Publishers
    //----------------------------------
    //mean_advertise = param_n.advertise<sensor_msgs::PointCloud2>("mean_map", 20);
    //var_advertise = param_n.advertise<sensor_msgs::PointCloud2>("var_map", 20);
    //----------------------------------
    // Services
    //----------------------------------
    //ros::ServiceServer service = param_n.advertiseService("suggestNextObservationLocation", suggestNextObservationLocation);

    module_init = false;
}


Cgmrf::~Cgmrf(){}


//-------------
// CALLBACKS
//-------------
void Cgmrf::sensorCallback(const olfaction_msgs::msg::GasSensor::SharedPtr msg)
{
    mutex_nose.lock();
    try
    {
        if(msg->raw_units == msg->UNITS_PPM)
        {
            //Normalize to [0,1]
            curr_reading = ((double)msg->raw - min_sensor_val)/ (max_sensor_val-min_sensor_val);
            if (curr_reading < 0)
                curr_reading = 0.0;
            if (curr_reading > 1)
                curr_reading = 1.0;
        }
        else if(msg->raw_units == msg->UNITS_OHM)
            curr_reading = (msg->raw_air - msg->raw)/msg->raw_air;  //range[0,1]
        else if(msg->raw_units == msg->UNITS_VOLT)
            curr_reading = msg->raw;
        else if(msg->raw_units == msg->UNITS_AMP)
            curr_reading = msg->raw;
    }catch(std::exception e){
        RCLCPP_ERROR(this->get_logger(), "[GMRF] Exception at new Obs: %s ", e.what() );
    }

    mutex_nose.unlock();


    //Get pose of the sensor in the /map reference system
    geometry_msgs::msg::TransformStamped t;
    std::string fromFrameRel = frame_id.c_str();
    std::string toFrameRel = msg->header.frame_id.c_str();
    try {
        t = tf_buffer_->lookupTransform(toFrameRel, fromFrameRel, tf2::TimePointZero);
    } catch (const tf2::TransformException & ex) {
        RCLCPP_ERROR(this->get_logger(), "Could not transform %s to %s: %s",
        toFrameRel.c_str(), fromFrameRel.c_str(), ex.what());
        return;
    }

    if (module_init)
    {
        //Current sensor pose
        float x_pos = t.transform.translation.x;
        float y_pos = t.transform.translation.y;

        //Add new observation to the map
        if (curr_reading <0 || curr_reading > 1)
        {
            RCLCPP_WARN(this->get_logger(), "[GMRF] Obs is out of bouns! %.2f [0,1]. Normalizing!", curr_reading);
            curr_reading = 1.0;
        }
        mutex_nose.lock();
        //ROS_WARN("[GMRF] New obs: %.2f at (%.2f,%.2f)", curr_reading,x_pos,y_pos);
        my_map->insertObservation_GMRF(curr_reading, x_pos, y_pos, GMRF_lambdaObs);
        mutex_nose.unlock();
    }
}



void Cgmrf::mapCallback(const nav_msgs::msg::OccupancyGrid::SharedPtr msg)
{
    RCLCPP_INFO(this->get_logger(), "%s - Got the map of the environment!", __FUNCTION__);
    occupancyMap = *msg;

    // Set GasMap dimensions as the OccupancyMap
    double map_min_x = msg->info.origin.position.x;
    double map_max_x = msg->info.origin.position.x + msg->info.width*msg->info.resolution;
    double map_min_y = msg->info.origin.position.y;
    double map_max_y = msg->info.origin.position.y + msg->info.height*msg->info.resolution;

    // Create GasMap
    //my_map = new CGMRF_map(occupancyMap, cell_size, GMRF_lambdaPrior,colormap,max_pclpoints_cell);
    my_map = new CGMRF_map(occupancyMap, cell_size, GMRF_lambdaPrior);
    RCLCPP_INFO(this->get_logger(), "[GMRF] GasGridMap initialized");

    module_init = true;
}


void Cgmrf::publishMaps()
{
    /*
    sensor_msgs::PointCloud2 meanPC, varPC;
    my_map->get_as_pointClouds(meanPC, varPC);
    meanPC.header.frame_id = frame_id.c_str();
    varPC.header.frame_id = frame_id.c_str();

    //ROS_INFO("[GMRF] Publishing maps!");
    mean_advertise.publish(meanPC);
    var_advertise.publish(varPC);
    */

    //SAVE AS CSV FILE
    if (output_csv_folder != "")
    {
        my_map->save_as_CSV(output_csv_folder);
    }
}


//----------------------------
//      MAIN
//----------------------------
int main(int argc, char **argv)
{
	rclcpp::init(argc, argv);

    // Greate a class object
    auto my_gmrf_map = std::make_shared<Cgmrf>();
	rclcpp::Rate loop_rate(my_gmrf_map->exec_freq);
	
	while (rclcpp::ok())
	{
		rclcpp::spin_some(my_gmrf_map);

		if (my_gmrf_map->module_init)
		{
			// Update and Publish maps
            my_gmrf_map->my_map->updateMapEstimation_GMRF(my_gmrf_map->GMRF_lambdaObsLoss);
            my_gmrf_map->publishMaps();
		}
        else
        {
			printf("[gmrf_node] Waiting for initialization (Map of environment).");
		}
        loop_rate.sleep();
	}
}