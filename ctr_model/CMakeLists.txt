cmake_minimum_required(VERSION 2.8.3)
project(ctr_model)

## Use C++14
add_definitions(--std=c++14)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  roscpp
  sensor_msgs
  tf
  tf2_ros
  urdf
  dynamic_reconfigure
  rospy
)

find_package(OMPL)

## Generate dynamic reconfigure parameters in the 'cfg' folder
generate_dynamic_reconfigure_options(
   cfg/MyParams.cfg
#   cfg/DynReconf2.cfg
)
###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS
    include
#  LIBRARIES
  CATKIN_DEPENDS 
    roscpp
    sensor_msgs
    tf
    tf2_ros
    urdf
#  DEPENDS 
)
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES anthropomorphic_arm_desrcription
#  CATKIN_DEPENDS urdf
#  DEPENDS system_lib
)
###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  /usr/include
  src/alglib
  include/alglib
)

## Declare a C++ executable
add_executable(${PROJECT_NAME}
  src/${PROJECT_NAME}_node.cpp
  src/CtrModel.cpp
  src/Calculation.cpp

  src/alglib/ap.cpp
  src/alglib/integration.cpp
  src/alglib/interpolation.cpp
  src/alglib/alglibinternal.cpp
  src/alglib/linalg.cpp
  src/alglib/alglibmisc.cpp
  src/alglib/specialfunctions.cpp
  src/alglib/solvers.cpp
  src/alglib/optimization.cpp 

)
add_dependencies(ctr_model ${PROJECT_NAME}_gencfg)
add_library(alglib
  src/alglib/alglibinternal.cpp
  src/alglib/alglibmisc.cpp
  src/alglib/ap.cpp
  src/alglib/linalg.cpp
  src/alglib/optimization.cpp
  src/alglib/solvers.cpp
)
target_include_directories(alglib PUBLIC include/alglib)
target_compile_definitions(alglib PUBLIC -DAE_CPU=AE_INTEL)
## Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
)


