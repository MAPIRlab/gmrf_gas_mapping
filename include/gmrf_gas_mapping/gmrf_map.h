
#include "rclcpp/rclcpp.hpp"
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <sensor_msgs/msg/point_cloud2.hpp>
#include <eigen3/Eigen/Sparse>
#include <fstream>      // std::ofstream

/*
#include <pcl/conversions.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
*/

#define NUM_CELL_TEMPLATES 200      //For plotting only
struct TRandomFieldCell
{
    double mean;
    double std;
};



/** GMRF class implementing the probability map and methods for insterting new observations and update the map */
class CGMRF_map: public rclcpp::Node
{
public:
    //CGMRF_map(const nav_msgs::msg::OccupancyGrid &oc_map, float cell_size, float m_lambdaPrior, std::string m_colormap, int max_points_cell);
    CGMRF_map(const nav_msgs::msg::OccupancyGrid &oc_map, float cell_size, float m_lambdaPrior);
    ~CGMRF_map();

    // insert new observation
    void insertObservation_GMRF(float normReading, float x_pos, float y_pos, float lambdaObs);

    // solves the minimum quadratic system to determine the new concentration of each cell
    void updateMapEstimation_GMRF(float lambdaObsLoss);

    // get mean and std maps as pointClouds (for visualization)
    // void get_as_pointClouds(sensor_msgs::msg::PointCloud2 &meanPC, sensor_msgs::msg::PointCloud2 &varPC);

    // save map as CSV file
    void save_as_CSV(std::string output_csv_file);

protected:
    std::vector<TRandomFieldCell>           m_map;
    float                                   m_x_min,m_x_max,m_y_min,m_y_max;
    float                                   m_resolution;       //cell_size
    size_t                                  m_size_x, m_size_y;
    size_t                                  N;                  //number of cells
    float                                   res_coef;

    float                                   GMRF_lambdaPrior;		//!< The information (Lambda) of fixed map constraints

    std::vector<Eigen::Triplet<double> >    H_prior;        // the prior part of H
    Eigen::VectorXd                         g;              // Gradient vector
    size_t                                  nPriorFactors;	// L
    size_t                                  nObsFactors;	// M
    size_t                                  nFactors;		// L+M
    std::multimap<size_t,size_t>            cell_interconnections;		//Store the interconnections (relations) of each cell with its neighbourds

    struct TobservationGMRF
    {
        float	obsValue;
        float	Lambda;
        bool	time_invariant;						//if the observation will lose weight (lambda) as time goes on (default false)
    };

    std::vector<std::vector<TobservationGMRF> > activeObs;		//Vector with the active observations and their respective Information


    /** Check if two cells of the gridmap (m_map) are connected, based on the provided occupancy gridmap*/
    bool exist_relation_between2cells(
        const nav_msgs::msg::OccupancyGrid *m_Ocgridmap,
        size_t cxo_min,
        size_t cxo_max,
        size_t cyo_min,
        size_t cyo_max,
        const size_t seed_cxo,
        const size_t seed_cyo,
        const size_t objective_cxo,
        const size_t objective_cyo
    );

    inline int x2idx(float x) const { return static_cast<int>( (x-m_x_min)/m_resolution ); }
    inline int y2idx(float y) const { return static_cast<int>( (y-m_y_min)/m_resolution ); }
    inline int xy2idx(float x,float y) const { return x2idx(x) + y2idx(y)*m_size_x; }

    // Visualization
    // pcl::PointCloud<pcl::PointXYZRGB> template_cells[NUM_CELL_TEMPLATES];
    // void init_pcl_templates(std::string colormap, int max_points_cell);
};
