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

#include "../mobile_robot.h"//机器人仿真模型文件

#include <iostream>
#include <fstream>

// Include file with properties
#include "../mobile_robot_wall_cts.h"//读取的一部分参数（eg. 机器人初始位置估计）

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
   * Linear system model      *
   ***************************/

  // Create the matrices A and B for the linear system model
  //x(k)=Ax(k-1)+Bu(k-1), x(k)为当前位置，u(k-1)速度量
  Matrix A(2,2);
  A(1,1) = 1.0;
  A(1,2) = 0.0;
  A(2,1) = 0.0;
  A(2,2) = 1.0;
  Matrix B(2,2);
  B(1,1) = cos(0.8);
  B(1,2) = 0.0;
  B(2,1) = sin(0.8);
  B(2,2) = 0.0;
//两个向量组合成一个vecter
  vector<Matrix> AB(2);
  AB[0] = A;
  AB[1] = B;

  // create gaussian
  //之后我们用均值向量μ以及一个协方差矩阵cov生成高斯分布以代表预测位置的不确定性。
  //均值由MU_SYSTEM_NOISE_X和MU_SYSTEM_NOISE_Y确定：
  ColumnVector sysNoise_Mu(2);
  sysNoise_Mu(1) = MU_SYSTEM_NOISE_X;
  sysNoise_Mu(2) = MU_SYSTEM_NOISE_Y;
  //将协方差矩阵选为对角矩阵，对角线上为MU_SYSTEM_NOISE_X和MU_SYSTEM_NOISE_Y：
  SymmetricMatrix sysNoise_Cov(2);
  sysNoise_Cov = 0.0;
  sysNoise_Cov(1,1) = SIGMA_SYSTEM_NOISE_X;
  sysNoise_Cov(1,2) = 0.0;
  sysNoise_Cov(2,1) = 0.0;
  sysNoise_Cov(2,2) = SIGMA_SYSTEM_NOISE_Y;
//均值与协方差矩阵一起定义二维高斯分布
  Gaussian system_Uncertainty(sysNoise_Mu, sysNoise_Cov);

  // create the model
  //创造一个线性条件概率密度函数（pdf）用以代表在给定当前机器人位置情况下，机器人预测位置的概率。
  //这个概率密度函数由两个参数构成：AB代表线性模型本身，高斯分布则代表系统模型的额外不确定度：
  LinearAnalyticConditionalGaussian sys_pdf(AB, system_Uncertainty);
  //终于我们通过系统的概率密度函数创造出整个预测系统的模型：
  LinearAnalyticSystemModelGaussianUncertainty sys_model(&sys_pdf);

  /*********************************
   * Initialise measurement model *
   ********************************/
//现在为了纠正机器人移动预测的位置，基于机器人距离墙的距离的测量值，建立测量更新模型．
//本例中使用线性测量更新模型，BFL包中已经包装好了。
//首先查看将测量值与机器人位置联系起来的数学公式：
//z(k+1)=Hx(k+1)
  // create matrix H for linear measurement model
  //建立矩阵Ｈ　
  Matrix H(1,2);
  double wall_ct = 2/(sqrt(pow(RICO_WALL,2.0) + 1));
  H = 0.0;
  H(1,1) = wall_ct * RICO_WALL;
  H(1,2) = 0 - wall_ct;
  //其中，wall_ct　是一个帮助常数，RICO_WALI是墙面的斜率．然后我们以传感器本身的
  //不确定性来建立一个高斯分布，这个高斯分布同样包含均值以及协方差矩阵cov
  // Construct the measurement noise (a scalar in this case)
  ColumnVector measNoise_Mu(1);
  measNoise_Mu(1) = MU_MEAS_NOISE;//噪声均值

  SymmetricMatrix measNoise_Cov(1);
  measNoise_Cov(1,1) = SIGMA_MEAS_NOISE;//噪声协方差矩阵
  Gaussian measurement_Uncertainty(measNoise_Mu, measNoise_Cov);
//现在可以创建测量可能的概率密度模型，这个概率密度模型由表示测量系统本身的矩阵Ｈ以及代表测量不确定性的
//高斯分布构成．最后根据这个概率密度模型我们可以构建测量更新模型：
  // create the model
  LinearAnalyticConditionalGaussian meas_pdf(H, measurement_Uncertainty);
  LinearAnalyticMeasurementModelGaussianUncertainty meas_model(&meas_pdf);


  /****************************
   * Linear prior DENSITY     先验分布*
   ***************************/
  //现在建立代表初始估计与不确定度的先验分布，同样包含了均值以及协方差如下所示：
   // Continuous Gaussian prior (for Kalman filters)
  ColumnVector prior_Mu(2);
  prior_Mu(1) = PRIOR_MU_X;
  prior_Mu(2) = PRIOR_MU_Y;
  SymmetricMatrix prior_Cov(2);
  prior_Cov(1,1) = PRIOR_COV_X;
  prior_Cov(1,2) = 0.0;
  prior_Cov(2,1) = 0.0;
  prior_Cov(2,2) = PRIOR_COV_Y;
  Gaussian prior(prior_Mu,prior_Cov);


  /******************************
   * Construction of the Filter构造滤波器 *
   ******************************/
  ExtendedKalmanFilter filter(&prior);


  /***************************
   * initialise MOBILE ROBOT *
   **************************/
  // Model of mobile robot in world with one wall
  // The model is used to simultate the distance measurements.
  //本文中的仿真机器人是在开放世界有一堵墙，机器人在这个开放世界中移动。
  //机器人的相关参数被定义在mobile_robot_wall_cts.h中（例如初始位置，仿真噪声等）。机器人仿真这样声明：
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
