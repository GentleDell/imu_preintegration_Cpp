#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "Lib_preintegration.h"
#include "se3.h"
#include "se3_ops.h"

#define _g 9.7964
#define pi 3.1415926

using namespace std;

typedef Matrix<double, 6, 1> Vector6d;

Eigen::Matrix3d I3x3 = Eigen::Matrix<double, 3, 3>::Identity();
Eigen::Matrix3d Z3x3 = Eigen::Matrix<double, 3, 3>::Identity();
Eigen::Vector3d I3(1.0, 1.0, 1.0);
Eigen::Vector3d Z3(0.0, 0.0, 0.0);

ImuPreintegration::ImuPreintegration()
{
    double default_gyro_noise = (2/180*pi);
    double default_acc_noise = 0.1;

    delta_tij = 0;
    delta_pij = Eigen::Vector3d(0,0,0);
    delta_vij = Eigen::Vector3d(0,0,0);
    bias_a = Eigen::Vector3d(0,0,0);
    bias_g = Eigen::Vector3d(0,0,0);
    bias_a_previous = Eigen::Vector3d(0,0,0);
    bias_g_previous = Eigen::Vector3d(0,0,0);
    delta_Rij = Eigen::Matrix<double, 3, 3>::Identity();
    ResiError << Z3,Z3,Z3;
    IMU_cov_ij = Eigen::Matrix<double, 9, 9>::Identity(); // simply initiate to identity

    // derivative init
    df_dx.DRij_Dbg = Z3x3;
    df_dx.Dpij_Dba = Z3x3;
    df_dx.Dpij_Dbg = Z3x3;
    df_dx.Dvij_Dba = Z3x3;
    df_dx.Dvij_Dbg = Z3x3;

    // covariance of imu measurement
    raw_acc_cov = default_acc_noise * default_acc_noise * Eigen::Matrix<double, 3, 3>::Identity();
    raw_gyro_cov = default_gyro_noise * default_gyro_noise * Eigen::Matrix<double, 3, 3>::Identity();
}

ImuPreintegration::ImuPreintegration(Eigen::Matrix3d acc_cov_input, Eigen::Matrix3d gyro_cov_input)
{
    delta_tij = 0;
    delta_pij = Eigen::Vector3d(0,0,0);
    delta_vij = Eigen::Vector3d(0,0,0);
    bias_a = Eigen::Vector3d(0,0,0);
    bias_g = Eigen::Vector3d(0,0,0);
    bias_a_previous = Eigen::Vector3d(0,0,0);
    bias_g_previous = Eigen::Vector3d(0,0,0);
    delta_Rij = Eigen::Matrix<double, 3, 3>::Identity();
    ResiError << Z3,Z3,Z3;
    IMU_cov_ij = Eigen::Matrix<double, 9, 9>::Identity(); // simply initiate to identity

    // derivative init
    df_dx.DRij_Dbg = Z3x3;
    df_dx.Dpij_Dba = Z3x3;
    df_dx.Dpij_Dbg = Z3x3;
    df_dx.Dvij_Dba = Z3x3;
    df_dx.Dvij_Dbg = Z3x3;

    // covariance of imu measurement
    raw_acc_cov = acc_cov_input;
    raw_gyro_cov = gyro_cov_input;
}

void ImuPreintegration::Preintegration( Eigen::Vector3d acc, Eigen::Vector3d gyro, double delta_t )
{
    delta_tij += delta_t;
    double t2_2 = delta_t*delta_t/2;
    SE3Quat exp_temp;

/// correct gyro measurement with bias
    Eigen::Vector3d gyroCorrect = gyro - bias_g;
    Eigen::Vector3d theta_gyro_corr = gyroCorrect*delta_t;
    Vector6d r_t;
    r_t << theta_gyro_corr, I3;
/// calculate the orientation change between consecutive time slots
    Eigen::Matrix3d Rk_k_1 = exp_temp.exp(r_t).rotation().matrix();
    Eigen::Matrix3d Rk_k_1_T = Rk_k_1.transpose();
/// calculate Jacobian —— Jr
    Eigen::Matrix3d Jrk_k_1 = Dexp( theta_gyro_corr ); // used in separating noise

/// correct acc measurement with bias
    Eigen::Vector3d accCorrect = acc - bias_a;
    Eigen::Matrix3d accCorrect_skew = skew(accCorrect);

/// incrementally update preintegration measurement noise covariance with new IMU measurements
/// supplementary fomula A.7, format is: [delta_hpi, delta_pose, delta_velocity]
    int size = 9;   // demension of Ak
    Eigen::MatrixXd Ak(size,size), Bk(size,3), Ck(size,3);

    Ak.topLeftCorner(size/3, size/3) = Rk_k_1_T;
    Ak.block(0,3,3,3) = Z3x3;  // from (0, 0)
    Ak.topRightCorner(size/3, size/3) = Z3x3;

    Ak.block(3,0,3,3) = -delta_Rij*accCorrect_skew*t2_2;   // there shuold be a 1/2?
    Ak.block(3,3,3,3) = I3x3;
    Ak.block(3,6,3,3) = I3x3*delta_t;

    Ak.bottomLeftCorner(size/3, size/3) = -delta_Rij*accCorrect_skew*delta_t;
    Ak.block(6,3,3,3) = Z3x3;  // from (0, 0)
    Ak.bottomRightCorner(size/3, size/3) = I3x3;

    Bk.topRows(3) = Z3x3;
    Bk.middleRows(3,3) = delta_Rij*t2_2;
    Bk.bottomRows(3) = delta_Rij*delta_t;

    Ck.topRows(3) = Jrk_k_1;
    Ck.middleRows(3,3) = Z3x3;
    Ck.bottomRows(3) = Z3x3;

/// SIGMA of pre-mea noise
    IMU_cov_ij = Ak*IMU_cov_ij*Ak.transpose() + Bk*raw_acc_cov*Bk.transpose() + Ck*raw_gyro_cov*Ck.transpose();  // covariance of IMU_ij residual errors

/// save previous derivative to bias of acc & gyro
    Eigen::Matrix3d DRij_Dbg = df_dx.DRij_Dbg;
    Eigen::Matrix3d Dpij_Dba = df_dx.Dpij_Dba;
    Eigen::Matrix3d Dpij_Dbg = df_dx.Dpij_Dbg;
    Eigen::Matrix3d Dvij_Dba = df_dx.Dvij_Dba;
    Eigen::Matrix3d Dvij_Dbg = df_dx.Dvij_Dbg;

/// incrementally update pre-mea derivative —— used for bias update
    df_dx.DRij_Dbg = Rk_k_1_T * DRij_Dbg - Jrk_k_1*delta_t;
    df_dx.Dpij_Dbg = Dpij_Dbg + Dvij_Dbg*delta_t - delta_Rij * accCorrect_skew*DRij_Dbg*t2_2;
    df_dx.Dpij_Dba = Dpij_Dba + Dvij_Dba * delta_t - delta_Rij*t2_2;
    df_dx.Dvij_Dbg = Dvij_Dbg - delta_Rij * accCorrect_skew*DRij_Dbg*delta_t;
    df_dx.Dvij_Dba = Dvij_Dba - delta_Rij * delta_t;

/// update preintegrated IMU measurements —— used to update NavState
    Eigen::Matrix3d R = delta_Rij;
    delta_Rij = delta_Rij*Rk_k_1;
//    o.phiv_ = SO3.log(o.R_);      // avoid the use of phiv
    delta_pij = delta_pij + delta_vij*delta_t + R*accCorrect*t2_2;    // missing 1/2 ?
    delta_vij = delta_vij + R*accCorrect*delta_t;

}

void ImuPreintegration::IMUResdual_Jaccobian(Matrix3d R_i, Matrix3d R_j, Vector3d p_i, Vector3d p_j, Vector3d v_i, Vector3d v_j)
{
    SE3Quat exp_temp;
    Vector6d r_t;
    Eigen::MatrixXd J_ResErro(9,24);   // "R v p" to "fai_i pi vi fai_j pj vj ba bg"
    Eigen::Vector3d g(0,0,-_g);
    Eigen::Vector3d delta_bg = bias_g_previous - bias_g;
    Eigen::Vector3d delta_ba = bias_a_previous - bias_a;
////  ResiError(10:15) = [delta_ba;delta_bg];     // donot consider the bias drift temporarily

    J_ResErro = Eigen::MatrixXd::Zero(9, 24);
    double t = delta_tij;
    double t22 = t*t/2;

/// incorporating Bias update to compensate the R V P
    Eigen::Vector3d DeltaR_from_bg = df_dx.DRij_Dbg*delta_bg;  // 3×1，so it's fai
    r_t << DeltaR_from_bg, I3;
    Eigen::Matrix3d correct_ResRij = delta_Rij*exp_temp.exp(r_t).rotation().matrix();
    Eigen::Vector3d correct_ResVij = delta_vij + df_dx.Dvij_Dba*delta_ba + df_dx.Dvij_Dbg*delta_bg;
    Eigen::Vector3d correct_Pij = delta_pij + df_dx.Dpij_Dba*delta_ba + df_dx.Dpij_Dbg*delta_bg;

/// residual error besed on formula (37) in the paper
/// format is  [R, v, p]
    Eigen::Quaterniond q_Rij_reserr(correct_ResRij.transpose() * R_i.transpose() * R_j);
    exp_temp.setRotation(q_Rij_reserr);
    ResiError.head(3) = exp_temp.log().head(3);     // translation will not affect the rotation
    ResiError.segment(3,3) = R_i.transpose() * (v_j - v_i - g*t) - correct_ResVij;
    ResiError.tail(3) = R_i.transpose() * (p_j - p_i - v_i*t - g*t22) - correct_Pij;  // g need 1/2 ?

/// obtain Jaccobian of residual error
/// state order is [fai, pose, velovcity]', derivative variances is [Fia_i, pi, vi, Fia_j, pj, vj, ba, bg]'
// Jaccobian of R to variant
    Eigen::Matrix3d Jr_resRij_inv = Dlog(ResiError.head(3));
    r_t << ResiError.head(3), I3;
    Eigen::Matrix3d RErrPhi = exp_temp.exp(r_t).rotation().matrix();
    Eigen::Matrix3d Jr_DeltaR_from_bg = Dexp(DeltaR_from_bg);
    J_ResErro.block(0,0,3,3) = -Jr_resRij_inv * R_j.transpose() * R_i;  // d ResErro_deltaRij, d fia_i
    J_ResErro.block(0,9,3,3) = Jr_resRij_inv;          // d ResErro_deltaRij, d fiaj
    J_ResErro.block(0,21,3,3) = -Jr_resRij_inv*RErrPhi.transpose()*Jr_DeltaR_from_bg*df_dx.DRij_Dbg;    // d ResErro_deltaRij, d bg

// Jaccobian of v to variant
    J_ResErro.block(3,0,3,3) = skew( R_i.transpose()*(v_j - v_i - g*t) );     // d ResErro_deltavij, d fai_i
    J_ResErro.block(3,6,3,3) = -R_i.transpose();   // d ResErro_deltavij, d vi
    J_ResErro.block(3,15,3,3) = R_i.transpose();   // d ResErro_deltavij, d vj
    J_ResErro.block(3,18,3,3) = -df_dx.Dvij_Dba;   // d ResErro_deltavij, d ba
    J_ResErro.block(3,21,3,3) = -df_dx.Dvij_Dbg;   // d ResErro_deltavij, d bg

// Jaccobian of p to variant
    J_ResErro.block(6,0,3,3) = skew( R_i.transpose()*(p_j - p_i - v_i*t - g*t22) );     // g shuold multiple 1/2
    J_ResErro.block(6,3,3,3) = -I3x3;
    J_ResErro.block(6,6,3,3) = -R_i.transpose() * t;
    //    J_ResErro.block(6,12,3,3) = R_i.transpose();
    J_ResErro.block(6,12,3,3) = R_i.transpose()*R_j;    // d ResErro_deltapij, d pj
    J_ResErro.block(6,18,3,3) = -df_dx.Dpij_Dba;     // d ResErro_deltapij, d ba
    J_ResErro.block(6,21,3,3) = -df_dx.Dpij_Dbg;     // d ResErro_deltapij, d bg
}

Eigen::VectorXd ImuPreintegration::preict(Eigen::Matrix3d R_i, Eigen::Vector3d p_i, Eigen::Vector3d v_i)
{
    double t22 = delta_tij*delta_tij/2;
    Eigen::Vector3d g(0,0,-_g);

    // transposition of formula (31) in the paper
    Eigen::Matrix3d R = R_i*delta_Rij;
    Eigen::Vector3d v = R_i*delta_vij + v_i + g*delta_tij;
    Eigen::Vector3d p = R_i*delta_pij + p_i + v_i*delta_tij + g*t22;    // missing 1/2?

    SE3Quat se3_Liegrp(R,p);
    Vector6d se3_Liealg = se3_Liegrp.log();
    Eigen::VectorXd state_j;
    state_j << se3_Liealg, v;   // R p v at j

    return state_j;
}

Eigen::Matrix3d ImuPreintegration::Dexp( Eigen::Vector3d theta )
{
    double phi = theta.norm();
    double phiinv = 1/phi;
    Eigen::Matrix3d phix, Jr;

    phix = skew(theta);

    if(phi < 0.001){
        Jr = I3x3 - 0.5*phix;
    }
    else{
        Jr = I3x3 - (1-cos(phi))*phiinv*phiinv*phix+(phi-sin(phi))*phiinv*phiinv*phiinv*phix*phix;
    }
    return Jr;
}

Eigen::Matrix3d ImuPreintegration::Dlog(Eigen::Vector3d theta )
{
    double phi = theta.norm();
    double phiinv = 1/phi;
    Eigen::Matrix3d phix,JrInv;

    phix = skew(theta);

    if(phi<1e-3){
        JrInv = I3x3 + 0.5*phix;
    }
    else{
        // error? according to last formula in supplementary
        // material, it shuold be + rather than - ?
        JrInv = I3x3 + 0.5*phix + (phiinv*phiinv-(1+cos(phi))*phiinv/(2*sin(phi)))*phix*phix;
    }
    return JrInv;
}
