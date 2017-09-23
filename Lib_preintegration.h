#ifndef LIBPREINTEGRATION_H_
#define LIBPREINTEGRATION_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;

class ImuPreintegration{

public:

    ImuPreintegration();    // default constructor

    ImuPreintegration( Eigen::Matrix3d acc_cov_input, Eigen::Matrix3d gyro_cov_input ); // constructor with cov input

    ~ImuPreintegration()
    {
        cout << "ImuPreintegration has been destructed!" << endl;
    }

    // preintegration & update covariance of Imu residual errors
    void Preintegration(Eigen::Vector3d acc, Eigen::Vector3d gyro, double delta_t );

    // feen state at i & j, obtain the Imu residual error and its Jaccobian to some variants
    void IMUResdual_Jaccobian(Eigen::Matrix3d R_i, Eigen::Matrix3d R_j, Eigen::Vector3d p_i, Eigen::Vector3d p_j,  Eigen::Vector3d v_i, Eigen::Vector3d v_j);

    // feed R p v in i, return predicted R p v in j (in the form of "se3 vel_vector")
    Eigen::VectorXd preict(Eigen::Matrix3d R_i, Eigen::Vector3d p_i, Eigen::Vector3d v_i);

    Eigen::Matrix3d Dexp( Eigen::Vector3d theta );
    Eigen::Matrix3d Dlog( Eigen::Vector3d theta );

    // after every optimiztion, we will have a new bias. To compensate residualerror essily, we need compute the delta_bias, so we record last bias
    void write_bias_previous()
    {
        bias_a_previous = bias_a;
        bias_g_previous = bias_g;
    }


private:

    double delta_tij;
    Eigen::Vector3d delta_pij, delta_vij, bias_a, bias_g, bias_a_previous, bias_g_previous;
    Eigen::Matrix3d delta_Rij;
    Eigen::Matrix3d raw_acc_cov, raw_gyro_cov;
    Eigen::MatrixXd IMU_cov_ij;     // SIGMA of pre-mea noise

    struct Derivative{
       Eigen::Matrix3d DRij_Dbg;
       Eigen::Matrix3d Dpij_Dba;
       Eigen::Matrix3d Dpij_Dbg;
       Eigen::Matrix3d Dvij_Dba;
       Eigen::Matrix3d Dvij_Dbg;
    } df_dx;

    Eigen::Matrix<double, 9, 1> ResiError;

};
#endif
