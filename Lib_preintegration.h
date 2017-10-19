#ifndef LIBPREINTEGRATION_H_
#define LIBPREINTEGRATION_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#define _g 1.1

using namespace std;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;

Eigen::Matrix3d I3x3 = Eigen::Matrix<double, 3, 3>::Identity();
Eigen::Matrix3d Z3x3 = Eigen::Matrix<double, 3, 3>::Zero();
Eigen::Vector3d I3(1.0, 1.0, 1.0);
Eigen::Vector3d Z3(0.0, 0.0, 0.0);


class ImuPreintegration{

public:

/// public para
    Eigen::Vector3d g;

/// public function
    ImuPreintegration();    // default constructor

    ImuPreintegration( Eigen::Matrix3d acc_cov_input, Eigen::Matrix3d gyro_cov_input ); // constructor with cov input

    ~ImuPreintegration()
    {
        cout << "ImuPreintegration has been destructed!" << endl;
    }

    ImuPreintegration( const ImuPreintegration& origin_IMUP );  // copy

    // reset the imupreintegraion to the original state
    void Reset();

    // preintegration & update covariance of Imu residual errors
    void Preintegration(Eigen::Vector3d acc, Eigen::Vector3d gyro, double delta_t );

    // feen state at i & j, obtain the Imu residual error and its Jaccobian to some variants
    void IMUResdual_Jaccobian(Eigen::Matrix3d R_i, Eigen::Matrix3d R_j, Eigen::Vector3d p_i, Eigen::Vector3d p_j,  Eigen::Vector3d v_i, Eigen::Vector3d v_j);

    // feed R p v in time i, return predicted R p v in time j (in the form of "se3 vel_vector")
    Vector9d predict(Eigen::Matrix3d R_i, Eigen::Vector3d p_i, Eigen::Vector3d v_i);

    // feed modifcaction, directly applied to this object.
    // If you want to generate a new object, just use copy construct function then call this func in new object.
    void update( Eigen::Matrix3d phiv_modify, Eigen::Vector3d p_modify, Eigen::Vector3d v_modify, Eigen::Vector3d bia_acc_modify, Eigen::Vector3d bia_gyro_modify );

    Eigen::Matrix3d Dexp( Eigen::Vector3d theta );
    Eigen::Matrix3d Dlog( Eigen::Vector3d theta );

    Eigen::Matrix3d Rij_read() { return delta_Rij; }
    Eigen::Vector3d pij_read() { return delta_pij; }
    Eigen::Vector3d vij_read() { return delta_vij; }
    double tij_read() { return delta_tij; }

    void Rij_write(Eigen::Matrix3d new_Rij) { delta_Rij = new_Rij; }
    void pij_write(Eigen::Vector3d new_pij) { delta_pij = new_pij; }
    void vij_write(Eigen::Vector3d new_vij) { delta_vij = new_vij; }

    // after every optimiztion, we will have a new bias. To compensate residualerror essily, we need compute the delta_bias, so we record last bias
    void write_bias_previous()
    {
        bias_a_previous = bias_a;
        bias_g_previous = bias_g;
    }

    inline Eigen::Quaterniond normalizeRotationQ(const Eigen::Quaterniond& quat)
    {
        Eigen::Quaterniond _quat(quat);
        if (_quat.w() < 0)
        {
            _quat.coeffs() *= -1;
        }
        return _quat.normalized();
    }

    inline Eigen::Matrix3d normalizeRotationM(const Eigen::Matrix3d& R)
    {
        Eigen::Quaterniond quat_R(R);
        return normalizeRotationQ(quat_R).toRotationMatrix();
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
