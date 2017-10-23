// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "types_IMU.h"

#include "../core/factory.h"
#include "../stuff/macros.h"

namespace g2o {

using namespace std;


VertexCamera::VertexCamera() : BaseVertex<9, Vector9d>() {
}

bool VertexCamera::read(std::istream& is)   // inputs are euler and translaiton & velocity
{
  Vector9d est;
  for (int i=0; i<9; i++)
    is  >> est[i];
  setEstimate(est);
  return true;
}

//bool VertexCamera::write(std::ostream& os) const {
////  SE3Quat cam2world(estimate().head(6).inverse());
////  for (int i=0; i<6; i++)
////    os << cam2world[i] << " ";
////  return os.good();
//    return false;
//}


EdgeIMUpreintegration::EdgeIMUpreintegration() : BaseBinaryEdge<9, Vector9d, VertexCamera, VertexCamera>() {
}

bool EdgeIMUpreintegration::read(std::istream& is){
  for (int i=0; i<=3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeIMUpreintegration::write(std::ostream& os) const {

  for (int i=0; i<=3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<=2; i++)
    for (int j=i; j<=2; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


/// obtain Jaccobian
/// state order is [fai, pose, velovcity]', derivative variances is [Fia_i, pi, vi, Fia_j, pj, vj, ba, bg]'
//// Jaccobian of R to variant
void EdgeIMUpreintegration::linearizeOplus()
{
    SE3Quat exp_temp;
    Vector6d r_t;
    double t = delta_tij;
    double t22 = t*t/2;
    // Eigen::MatrixXd J_ResErro(9,24);   // "R v p" to "fai_i pi vi fai_j pj vj ba bg"
    Eigen::Vector3d delta_bg = bias_g_previous - bias_g;
    Eigen::Vector3d DeltaR_from_bg = df_dx.DRij_Dbg*delta_bg;  // 3×1，so it's fai

    VertexCamera * vj = static_cast<VertexCamera *>(_vertices[1]);
    Vector9d state_j = vj->estimate();

    VertexCamera * vi = static_cast<VertexCamera *>(_vertices[0]);
    Vector9d state_i = vi->estimate();

    Quaterniond Q;
    Q = AngleAxisd(state_i[0], Vector3d::UnitX()) * AngleAxisd(state_i[1], Vector3d::UnitY()) * AngleAxisd(state_i[2], Vector3d::UnitZ());
    Matrix3d R_i = Q.toRotationMatrix();
    Vector3d p_i = state_i.segment(3,3);
    Vector3d v_i = state_i.tail(3);

    Q = AngleAxisd(state_j[0], Vector3d::UnitX()) * AngleAxisd(state_j[1], Vector3d::UnitY()) * AngleAxisd(state_j[2], Vector3d::UnitZ());
    Matrix3d R_j = Q.toRotationMatrix();
    Vector3d p_j = state_j.segment(3,3);
    Vector3d v_j = state_j.tail(3);

    // Jaccobian of R to variant
    Eigen::Matrix3d Jr_resRij_inv = Dlog(ResiError.head(3));
    r_t << ResiError.head(3), I3;
    Eigen::Matrix3d RErrPhi = exp_temp.exp(r_t).rotation().matrix();
    Eigen::Matrix3d Jr_DeltaR_from_bg = Dexp(DeltaR_from_bg);
    _jacobianOplusXi.block(0,0,3,3) = -Jr_resRij_inv * R_j.transpose() * R_i;  // d ResErro_deltaRij, d fia_i
    _jacobianOplusXi.block(0,9,3,3) = Jr_resRij_inv;          // d ResErro_deltaRij, d fiaj
    _jacobianOplusXi.block(0,21,3,3) = -Jr_resRij_inv*RErrPhi.transpose()*Jr_DeltaR_from_bg*df_dx.DRij_Dbg;    // d ResErro_deltaRij, d bg

    // Jaccobian of v to variant
    _jacobianOplusXi.block(3,0,3,3) = skew( R_i.transpose()*(v_j - v_i - g*t) );     // d ResErro_deltavij, d fai_i
    _jacobianOplusXi.block(3,6,3,3) = -R_i.transpose();   // d ResErro_deltavij, d vi
    _jacobianOplusXi.block(3,15,3,3) = R_i.transpose();   // d ResErro_deltavij, d vj
    _jacobianOplusXi.block(3,18,3,3) = -df_dx.Dvij_Dba;   // d ResErro_deltavij, d ba
    _jacobianOplusXi.block(3,21,3,3) = -df_dx.Dvij_Dbg;   // d ResErro_deltavij, d bg

    // Jaccobian of p to variant
    _jacobianOplusXi.block(6,0,3,3) = skew( R_i.transpose()*(p_j - p_i - v_i*t - 0.5*g*t22) );
    _jacobianOplusXi.block(6,3,3,3) = -I3x3;
    _jacobianOplusXi.block(6,6,3,3) = -R_i.transpose() * t;
    // J_ResErro.block(6,12,3,3) = R_i.transpose();
    _jacobianOplusXi.block(6,12,3,3) = R_i.transpose()*R_j;    // d ResErro_deltapij, d pj
    _jacobianOplusXi.block(6,18,3,3) = -df_dx.Dpij_Dba;     // d ResErro_deltapij, d ba
    _jacobianOplusXi.block(6,21,3,3) = -df_dx.Dpij_Dbg;     // d ResErro_deltapij, d bg


    _jacobianOplusXi = _jacobianOplusXj;
}

void EdgeIMUpreintegration::IMUResdual(Vector9d state_i, Vector9d state_j)
{
    Quaterniond Q;
    Q = AngleAxisd(state_i[0], Vector3d::UnitX()) * AngleAxisd(state_i[1], Vector3d::UnitY()) * AngleAxisd(state_i[2], Vector3d::UnitZ());
    Matrix3d R_i = Q.toRotationMatrix();
    Vector3d p_i = state_i.segment(3,3);
    Vector3d v_i = state_i.tail(3);

    Q = AngleAxisd(state_j[0], Vector3d::UnitX()) * AngleAxisd(state_j[1], Vector3d::UnitY()) * AngleAxisd(state_j[2], Vector3d::UnitZ());
    Matrix3d R_j = Q.toRotationMatrix();
    Vector3d p_j = state_j.segment(3,3);
    Vector3d v_j = state_j.tail(3);

    SE3Quat exp_temp;
    Vector6d r_t;
    Eigen::MatrixXd J_ResErro(9,24);   // "R v p" to "fai_i pi vi fai_j pj vj ba bg"
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
    ResiError.tail(3) = R_i.transpose() * (p_j - p_i - v_i*t - 0.5*g*t22) - correct_Pij;
}

void EdgeIMUpreintegration::Data_input(Eigen::Matrix3d deltaRij, Eigen::Vector3d deltapij, Eigen::Vector3d deltavij, Eigen::Vector3d biasa, Eigen::Vector3d biasg, Eigen::Vector3d biasa_previous, Eigen::Vector3d biasg_previous, double deltatij)
{
    delta_Rij = deltaRij;
    delta_pij = deltapij;
    delta_vij = deltavij;
    delta_tij = deltatij;
    bias_a = biasa;
    bias_g = biasg;
    bias_a_previous = biasa_previous;
    bias_g_previous = biasg_previous;
}

void EdgeIMUpreintegration::df_dx_input(Eigen::Matrix3d DRij_Dbg, Eigen::Matrix3d Dpij_Dba, Eigen::Matrix3d Dpij_Dbg, Eigen::Matrix3d Dvij_Dba, Eigen::Matrix3d Dvij_Dbg)
{
    df_dx.DRij_Dbg = DRij_Dbg;
    df_dx.Dpij_Dba = Dpij_Dba;
    df_dx.Dpij_Dbg = Dpij_Dbg;
    df_dx.Dvij_Dba = Dvij_Dba;
    df_dx.Dvij_Dbg = Dvij_Dbg;
}

Eigen::Matrix3d EdgeIMUpreintegration::Dexp( Eigen::Vector3d theta )
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


Eigen::Matrix3d EdgeIMUpreintegration::Dlog(Eigen::Vector3d theta )
{
    double phi = theta.norm();
    double phiinv = 1/phi;
    Eigen::Matrix3d phix,JrInv;

    phix = skew(theta);

    if(phi<1e-3){
        JrInv = I3x3 + 0.5*phix;
    }
    else{
        JrInv = I3x3 + 0.5*phix + (phiinv*phiinv-(1+cos(phi))*phiinv/(2*sin(phi)))*phix*phix;
    }
    return JrInv;
}


} // end namespace
