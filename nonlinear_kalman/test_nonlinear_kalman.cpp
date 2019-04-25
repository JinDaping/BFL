// $Id: kalman_mobile.cpp 5925 2006-03-14 21:23:49Z tdelaet $
// Copyright (C) 2006 Tinne De Laet <first dot last at mech dot kuleuven dot be>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.

/* Demonstration program for the Bayesian Filtering Library.
   Mobile robot localization with respect to wall with different possibilities for filter
*/


#include <filter/extendedkalmanfilter.h>

#include <model/linearanalyticsystemmodel_gaussianuncertainty.h>
#include <model/linearanalyticmeasurementmodel_gaussianuncertainty.h>

#include <pdf/analyticconditionalgaussian.h>
#include <pdf/linearanalyticconditionalgaussian.h>
#include "../nonlinearanalyticconditionalgaussianmobile.h"//added

#include "../mobile_robot.h"

#include <iostream>
#include <fstream>

// Include file with properties
#include "../mobile_robot_wall_cts.h"

using namespace MatrixWrapper;
using namespace BFL;
using namespace std;



/* The purpose of this program is to construct a kalman filter for the problem
  of localisation of a mobile robot equipped with an ultrasonic sensor.
  In this case the orientation is known, which simplifies the model considerably:
  The system model will become linear.
  The ultrasonic measures the distance to the wall (it can be switched off:
  see mobile_robot_wall_cts.h)

  The necessary SYSTEM MODEL is:

  x_k      = x_{k-1} + v_{k-1} * cos(theta) * delta_t
  y_k      = y_{k-1} + v_{k-1} * sin(theta) * delta_t

  The used MEASUREMENT MODEL:
  measuring the (perpendicular) distance z to the wall y = ax + b

  set WALL_CT = 1/sqrt(pow(a,2) + 1)
  z = WALL_CT * a * x - WALL_CT * y + WALL_CT * b + GAUSSIAN_NOISE
  or Z = H * X_k + J * U_k

  where

  H = [ WALL_CT * a       - WALL_CT      0 ]
  and GAUSSIAN_NOISE = N((WALL_CT * b), SIGMA_MEAS_NOISE)

*/


int main(int argc, char** argv)
{
  cerr << "==================================================" << endl
       << "Test of kalman filter" << endl
       << "Mobile robot localisation example" << endl
       << "==================================================" << endl;


  /****************************
   * NonLinear system model      *
   ***************************/
//均值以及协方差矩阵的值均定义在mobie_robot_wall_cts.h
//
  // create gaussian
  //预测噪声均值mu
  ColumnVector sys_noise_Mu(STATE_SIZE);
  sys_noise_Mu(1) = MU_SYSTEM_NOISE_X;
  sys_noise_Mu(2) = MU_SYSTEM_NOISE_Y;
  sys_noise_Mu(3) = MU_SYSTEM_NOISE_THETA;
//预测噪声协方差矩阵cov,为对角矩阵
  SymmetricMatrix sys_noise_Cov(STATE_SIZE);
  sys_noise_Cov = 0.0;
  sys_noise_Cov(1,1) = SIGMA_SYSTEM_NOISE_X;
  sys_noise_Cov(1,2) = 0.0;
  sys_noise_Cov(1,3) = 0.0;
  sys_noise_Cov(2,1) = 0.0;
  sys_noise_Cov(2,2) = SIGMA_SYSTEM_NOISE_Y;
  sys_noise_Cov(2,3) = 0.0;
  sys_noise_Cov(3,1) = 0.0;
  sys_noise_Cov(3,2) = 0.0;
  sys_noise_Cov(3,3) = SIGMA_SYSTEM_NOISE_THETA;
//构建一个三维高斯分布以代表预测模型的不确定
  Gaussian system_Uncertainty(sys_noise_Mu, sys_noise_Cov);

  // create the model
  //非线性概率密度函数//高斯不确定度的解析概率密度模型
  NonLinearAnalyticConditionalGaussianMobile sys_pdf(system_Uncertainty);
  AnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);

  /*********************************
   * Initialise measurement model 建立测量更新模型*
   ********************************/

  // create matrix H for linear measurement model
  //因为状态向量多了一个机器人方向，所以测量更新模型多了额外的０，其他步骤与之前的例子一样
  double wall_ct = 2/(sqrt(pow(RICO_WALL,2.0) + 1));
  Matrix H(MEAS_SIZE,STATE_SIZE);
  H = 0.0;
  H(1,1) = wall_ct * RICO_WALL;
  H(1,2) = 0 - wall_ct;
  H(1,3) = 0.0;
  // Construct the measurement noise (a scalar in this case)
  ColumnVector meas_noise_Mu(MEAS_SIZE);
  meas_noise_Mu(1) = MU_MEAS_NOISE;

  SymmetricMatrix meas_noise_Cov();
  meas_noise_Cov(1,1) = SIGMA_MEAS_NOISE;
  Gaussian measurement_Uncertainty(meas_noise_Mu, meas_noise_Cov);

  // create the measurement model
  LinearAnalyticConditionalGaussian meas_pdf(H, measurement_Uncertainty);
  LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);

  /****************************
   * Linear prior DENSITY    先验分布 *
   ***************************/
   // Continuous Gaussian prior (for Kalman filters)
   //与上一个相比，鲜艳分布增加一个方向变量：
  ColumnVector prior_Mu(STATE_SIZE);
  prior_Mu(1) = PRIOR_MU_X;
  prior_Mu(2) = PRIOR_MU_Y;
  prior_Mu(3) = PRIOR_MU_THETA;
  //协方差矩阵是对角矩阵
  SymmetricMatrix prior_Cov(STATE_SIZE);
  prior_Cov(1,1) = PRIOR_COV_X;
  prior_Cov(1,2) = 0.0;
  prior_Cov(1,3) = 0.0;
  prior_Cov(2,1) = 0.0;
  prior_Cov(2,2) = PRIOR_COV_Y;
  prior_Cov(2,3) = 0.0;
  prior_Cov(3,1) = 0.0;
  prior_Cov(3,2) = 0.0;
  prior_Cov(3,3) = PRIOR_COV_THETA;
  //均值与协方差矩阵构建成一个三维的高斯分布
  Gaussian prior_cont(prior_Mu,prior_Cov);

  /******************************
   * Construction of the Filter *
   ******************************/
  //扩展卡尔曼滤波器估计移动机器人的位置
  ExtendedKalmanFilter filter(&prior_cont);

  /***************************
   * initialise MOBILE ROBOT *
   **************************/
  // Model of mobile robot in world with one wall
  // The model is used to simultate the distance measurements.
  MobileRobot mobile_robot;
  ColumnVector input(2);
  input(1) = 0.1;
  input(2) = 0.0;




  /*******************
   * ESTIMATION LOOP *
   *******************/
  cout << "MAIN: Starting estimation" << endl;
  unsigned int time_step;
  for (time_step = 0; time_step < NUM_TIME_STEPS-1; time_step++)
    {
      // DO ONE STEP WITH MOBILE ROBOT
      mobile_robot.Move(input);

      // DO ONE MEASUREMENT
      ColumnVector measurement = mobile_robot.Measure();

      // UPDATE FILTER
      filter.Update(&sys_model,input,&meas_model,measurement);
      //filter.Update(&sys_model,input);

    } // estimation loop



  Pdf<ColumnVector> * posterior = filter.PostGet();
  cout << "After " << time_step+1 << " timesteps " << endl;
  cout << " Posterior Mean = " << endl << posterior->ExpectedValueGet() << endl
       << " Covariance = " << endl << posterior->CovarianceGet() << "" << endl;


  cout << "======================================================" << endl
       << "End of the Kalman filter for mobile robot localisation" << endl
       << "======================================================"
       << endl;


  return 0;
}
