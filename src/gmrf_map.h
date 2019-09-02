
#include "ros/ros.h"
#include <nav_msgs/OccupancyGrid.h>
#include <sensor_msgs/PointCloud2.h>
#include <Eigen/Sparse>
#include <fstream>      // std::ofstream
#include <pcl/conversions.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include "std_msgs/Float32MultiArray.h"
#include "std_msgs/MultiArrayDimension.h"
#include <vector>
#include <array>



#define NUM_CELL_TEMPLATES 200      //For plotting only
struct TRandomFieldCell
{
    double mean;
    double std;
};



/** GMRF class implementing the probability map and methods for insterting new observations and update the map */
class CGMRF_map
{
public:
    CGMRF_map(const nav_msgs::OccupancyGrid &oc_map, float cell_size, float m_lambdaPrior, std::string m_colormap, int max_points_cell);
    ~CGMRF_map();
    //insert new observation
    void  insertObservation_GMRF(float normReading, float x_pos, float y_pos, float lambdaObs);
    //solves the minimum quadratic system to determine the new concentration of each cell
    void  updateMapEstimation_GMRF(float lambdaObsLoss);
    //GET mean and std maps
    void get_as_pointClouds(sensor_msgs::PointCloud2 &meanPC, sensor_msgs::PointCloud2 &varPC);
    //Stores map as CSV to specified file
    void store_as_CSV(std::string output_csv_file);
    //Writes raw mean and std values to matrix
    void write_raw_values(std_msgs::Float32MultiArray &meanRAW, std_msgs::Float32MultiArray &varRAW);

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
        const nav_msgs::OccupancyGrid *m_Ocgridmap,
        size_t cxo_min,
        size_t cxo_max,
        size_t cyo_min,
        size_t cyo_max,
        const size_t seed_cxo,
        const size_t seed_cyo,
        const size_t objective_cxo,
        const size_t objective_cyo
    );

    inline int   x2idx(float x) const { return static_cast<int>( (x-m_x_min)/m_resolution ); }
    inline int   y2idx(float y) const { return static_cast<int>( (y-m_y_min)/m_resolution ); }
    inline int   xy2idx(float x,float y) const { return x2idx(x) + y2idx(y)*m_size_x; }

    //Visualization
    pcl::PointCloud<pcl::PointXYZRGB> template_cells[NUM_CELL_TEMPLATES];
    void init_pcl_templates(std::string colormap, int max_points_cell);

};
