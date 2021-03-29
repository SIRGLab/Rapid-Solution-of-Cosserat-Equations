# Rapid-Solution-of-Cosserat-Equations

Repo for rapid modeling of a Cosserat Rod in Matlab. It implements the methodology presented in the following publication:
M. Khadem, L. Da Cruz and C. Bergeles, "Rapid Solution of Cosserat Rod Equations via a Nonlinear Partial Observer," 2021 IEEE International Conference on Robotics and Automation  (ICRA 2021).

This repo uses C++ in ROS environment to rapidly simulate motion a Cosserat rod. You need to copy the 'ctr_model' folder into the src folder in your workspace. Then you can run it by the following command:  roslaunch ctr_model ctr_model.launch

If you enjoy this repository and use it, please cite our paper

```
@inproceedings{balintICRA21,
  title={Rapid Solution of Cosserat Rod Equations via a Nonlinear Partial Observer},
  author={Thamo, Balint and Dhalival, Kev and Khadem, Mohsen},
  booktitle={2021 IEEE International Conference on Robotics and Automation  (ICRA 2021)},
  pages={},
  year={2021},
  organization={IEEE}
}
```
