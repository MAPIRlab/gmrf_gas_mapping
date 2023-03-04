//----------------------------------------------------------------
// Simple ROS2 publisher of olfaction measurements from CSV file.
// Used as a simple test of GMRF_node
//-----------------------------------------------------
#include <rclcpp/rclcpp.hpp>
#include "std_msgs/msg/float32.hpp"
#include "olfaction_msgs/msg/observation.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

class Gmrf_pub: public rclcpp::Node
{
public:
    Gmrf_pub();
    ~Gmrf_pub();
    double exec_freq;
    bool publish_measurement();

protected:
    //Subscriptions & Publishers
    rclcpp::Publisher<olfaction_msgs::msg::Observation>::SharedPtr pub_;

    // params
    std::string frame_id, topic, input_csv_file;

    // CSV file 
    std::fstream file;
    std::string line, word;
    std::vector<std::string> row;
};