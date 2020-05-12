/* Kinematics.h
  
  Calculations of kinematics to convert from joint angles and coordinates for robot arm.
  This code is created by TheFLash/YWW.
*/

#ifndef Kinematics_h
#define Kinematics_h


  
 void InverseK(float* Xik, float* Jik);
 void ForwardK(float* Jfk, float* Xfk);
    

#endif
