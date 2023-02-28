import os

from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, SetEnvironmentVariable
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory


def generate_launch_description():
    # Get the launch directory
    my_dir = get_package_share_directory('gmrf_gas_mapping')

    # common variables
    use_sim_time = True
    remappings=[]
    params_yaml_file = os.path.join(my_dir, 'launch', 'gmrf_params.yaml')
    map_file = os.path.join(my_dir, 'maps', 'empty_map.yaml')
    
    logger = LaunchConfiguration("log_level")
    
    return LaunchDescription([
        # Set env var to print messages to stdout immediately
        SetEnvironmentVariable('RCUTILS_LOGGING_BUFFERED_STREAM', '1'),
        
        DeclareLaunchArgument(
            "log_level",
            default_value=["info"],  #debug, info
            description="Logging level",
            ),

        # MAP_SERVER
        Node(
            package='nav2_map_server',
            executable='map_server',
            name='map_server',
            output='screen',
            parameters=[{'use_sim_time': use_sim_time},
                        {'yaml_filename' : map_file}],
            remappings=remappings
            ),
        
        # GMRF
        Node(
            package='gmrf_gas_mapping',
            executable='gmrf_node',
            name='gmrf_node',
            output='screen',
            parameters=[params_yaml_file],
            remappings=remappings
            ),

        # LIFECYCLE MANAGER
        Node(
            package='nav2_lifecycle_manager',
            executable='lifecycle_manager',
            name='lifecycle_manager_localization',
            output='screen',
            parameters=[{'use_sim_time': use_sim_time},
                        {'autostart': True},
                        {'node_names': ['map_server']}
                       ]
            )
    ])
