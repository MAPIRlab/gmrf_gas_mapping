#include "gmrf_gas_mapping/gmrf_simple_publisher.hpp"


using std::placeholders::_1;

// -------------
// Gmrf_pub
//--------------
Gmrf_pub::Gmrf_pub(): Node("GMRF_pub")
{    
    RCLCPP_INFO(this->get_logger(), "[GMRF_pub] Loading Node parameters");

    //-----------------------------
    // Declare and Load Parameters
    //-----------------------------    
    frame_id = this->declare_parameter<std::string>("frame_id", "map");
    this->get_parameter("frame_id", frame_id);

    topic = this->declare_parameter<std::string>("topic", "/noTopic");
    this->get_parameter("topic", topic);

    input_csv_file = this->declare_parameter<std::string>("input_csv_file", "/noFile");
    this->get_parameter("input_csv_file", input_csv_file);

    exec_freq = this->declare_parameter<double>("exec_freq", 1.0);
    this->get_parameter("exec_freq", exec_freq);

    //----------------------------------
    // Publisher
    //----------------------------------
    pub_ = this->create_publisher<olfaction_msgs::msg::Observation>(topic, 10);

    // Open file
    file.open(input_csv_file, std::ios::in);
}

Gmrf_pub::~Gmrf_pub(){
    file.close();
}

bool Gmrf_pub::publish_measurement()
{
    // Read new line ( X Y ppm)
    if(file.is_open())
	{
		if(std::getline(file, line))
		{
			row.clear();
 
			std::stringstream my_str(line);
 
			while(std::getline(my_str, word, ','))
				row.push_back(word);
		}
        else
        {
            RCLCPP_INFO(this->get_logger(), "[GMRF_pub] End Of File: %s", input_csv_file.c_str());
            return false;
        }
	}
	else{
        RCLCPP_INFO(this->get_logger(), "[GMRF_pub] Error opening CSV file: %s", input_csv_file.c_str());
        return false;
    }
		
    // Compose msg
    olfaction_msgs::msg::Observation msg;
    msg.header.stamp = this->now();
    msg.header.frame_id = frame_id;

    msg.position.x = std::atof(row[0].c_str())*0.3;
    msg.position.y = std::atof(row[1].c_str())*0.3;
    msg.type = 1;   // gas only
    msg.gas = std::atof(row[2].c_str());
    
    // publish
    RCLCPP_INFO(this->get_logger(), "[GMRF_pub] Publishing new Obs [x,y,ppm]=[%.2f, %.2f, %.2f]", msg.position.x, msg.position.y, msg.gas);
    pub_->publish(msg);
    return true;
}



// MAIN
int main(int argc, char * argv[])
{
    sleep(5);

    rclcpp::init(argc, argv);  

    //create object
    std::shared_ptr<Gmrf_pub> mySP = std::make_shared<Gmrf_pub>();
    
    rclcpp::Rate loop_rate(mySP->exec_freq);
    bool file_reading = true;
    while (rclcpp::ok() && file_reading)
    {
        rclcpp::spin_some(mySP);
        file_reading = mySP->publish_measurement();
        loop_rate.sleep();
    }

    rclcpp::shutdown();
    return 0;
}