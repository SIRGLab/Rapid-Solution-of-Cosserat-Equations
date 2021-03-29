#include <ros/ros.h>
#include "ctr_model/CtrModel.hpp"

int main(int argc, char** argv)
{
  ros::init(argc, argv, "ctr_model");
  ros::NodeHandle nodeHandle("~");

  ctr_model::CtrModel CtrModel(nodeHandle);

  ros::spin();
  return 0;
}
