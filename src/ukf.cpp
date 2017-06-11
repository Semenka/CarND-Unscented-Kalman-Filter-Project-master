#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_ << 0, 0, 0, 0, 0;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  L = MatrixXd(2,2);
  L <<
  std_laspx_*std_laspx_,0,
  0,std_laspy_*std_laspy_;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  R = MatrixXd(3,3);
  R <<
  std_radr_*std_radr_,0,0,
  0,std_radphi_*std_radphi_,0,
  0,0,std_radrd_*std_radrd_;

  ///* State dimension
  n_x_=5;

  ///* Augmented state dimension
  n_aug_=7;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  ///* Sigma point spreading parameter
  lambda_=3 - n_aug_;

  ///* Weights of sigma points
  weights_=VectorXd(2*n_aug_+1);


  is_initialized_=false;



  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "Radar initialization"<< endl;
      float rho = meas_package.raw_measurements_[0]; // range
      float phi = meas_package.raw_measurements_[1]; // bearing
      float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho
      // Coordinates convertion from polar to cartesian
      float x = rho * cos(phi);
      float y = rho * sin(phi);
      x_ << x, y, rho_dot,0,0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "Laser initialization"<< endl;
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0,0;

    }


    // done initializing, no need to predict or update
    weights_(0)=lambda_/(lambda_+n_aug_);
    for (int i=1;i<2 * n_aug_ + 1;i++)
    {
      weights_(i)=1/(2*(lambda_+n_aug_));
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar updates

    UpdateRadar(meas_package);
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {

    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */


  //1.Generate sigma points
  //create sigma point matrix
  //MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  //MatrixXd A = P_.llt().matrixL();
  //Xsig.col(0)  = x_;

  //set remaining sigma points
  //for (int i = 0; i < n_x_; i++)
  //{
  //  Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
  //  Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  //}
  //2.Augmentation
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //3.Sigma points prediction


  VectorXd x_pred = VectorXd(n_x_);

  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
    x_aug=Xsig_aug.col(i);
    double v=x_aug(2);
    double phi=x_aug(3);
    double phi_dot=x_aug(4);
    double v_a=x_aug(5);
    double v_phi=x_aug(6);
    if (phi_dot!=0){
            x_pred(0)=x_aug(0)+v*(sin(phi+phi_dot*delta_t)-sin(phi))/phi_dot+delta_t*delta_t*cos(phi)*v_a/2;
            x_pred(1)=x_aug(1)+v*(-cos(phi+phi_dot*delta_t)+cos(phi))/phi_dot+delta_t*delta_t*sin(phi)*v_a/2;
            x_pred(2)=x_aug(2)+delta_t*v_a;
            x_pred(3)=x_aug(3)+phi_dot*delta_t+delta_t*delta_t*v_phi/2;
            x_pred(4)=x_aug(4)+delta_t*v_phi;

    }
    else {
            x_pred(0)=x_aug(0)+v*cos(phi)*delta_t+delta_t*delta_t*cos(phi)*v_a/2;
            x_pred(1)=x_aug(1)+v*sin(phi)*delta_t+delta_t*delta_t*cos(phi)*v_a/2;
            x_pred(2)=x_aug(2)+delta_t*v_a;
            x_pred(3)=x_aug(3)+delta_t*delta_t*v_phi/2;
            x_pred(4)=x_aug(4)+delta_t*v_phi;
    }

    Xsig_pred_.col(i) = x_pred;
  }


  //Calculation of mean and variance

  //predict state mean
  x_ = Xsig_pred_ * weights_;
  P_.fill(0.0);
  //predict state covariance matrix

  

  
  for (int i=0;i<2 * n_aug_ + 1;i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
      while (x_diff(1)<-M_PI) x_diff(1)+=2.*M_PI;

      P_=P_+(weights_(i)*x_diff*x_diff.transpose());
  }


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
      //transform sigma points into measurement space
  int n_z = 2;
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = Zsig*weights_;

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  for (int i=0;i<2 * n_aug_ + 1;i++)
  {
    S=S+weights_(i)*(Zsig.col(i)-z_pred)*(Zsig.col(i)-z_pred).transpose();
  }
  S=S+L;

  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
   //calculate cross correlation matrix
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
    Tc=Tc+weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
  }//end of for loop

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  x_=x_+K*(z-z_pred);
  P_=P_-K*S*K.transpose();
  NIS_lidar_=(z-z_pred).transpose()*S.inverse()*(z-z_pred);
  cout << "NIS_lidar_" << NIS_lidar_<< endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    //transform sigma points into measurement space
  int n_z = 3;
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);


  for (int i=0;i<2 * n_aug_ + 1;i++)
  {
      double px=Xsig_pred_.col(i)(0);
      double py=Xsig_pred_.col(i)(1);
      double v=Xsig_pred_.col(i)(2);
      double phi=Xsig_pred_.col(i)(3);
      double phi_dot=Xsig_pred_.col(i)(4);

      const double  eps=0.001;
      Zsig.col(i)(0)=sqrt(px*px+py*py);
      Zsig.col(i)(1)=atan2(py,px);
      Zsig.col(i)(2)=(px*cos(phi)*v+py*sin(phi)*v)/std::max(eps, Zsig.col(i)(0));

      z_pred=z_pred+weights_(i)*Zsig.col(i);

  }
  for (int i=0;i<2 * n_aug_ + 1;i++){
        S=S+weights_(i)*(Zsig.col(i)-z_pred)*(Zsig.col(i)-z_pred).transpose();
  }

  S=S+R;

  //create example vector for incoming radar measurement
  //set measurement dimension, radar can measure r, phi, and r_dot

  VectorXd z = meas_package.raw_measurements_;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
   //calculate cross correlation matrix
  for(int i=0;i<2 * n_aug_ + 1;i++)
  {
    Tc=Tc+weights_(i)*(Xsig_pred_.col(i)-x_)*(Zsig.col(i)-z_pred).transpose();
  }//end of for loop

  //calculate Kalman gain K;
  MatrixXd K =Tc*S.inverse();

    //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_=x_+K*z_diff;
  P_=P_-K*S*K.transpose();
  NIS_radar_=z_diff.transpose()*S.inverse()*z_diff;
  cout << "NIS_radar_" << NIS_radar_<< endl;

}
