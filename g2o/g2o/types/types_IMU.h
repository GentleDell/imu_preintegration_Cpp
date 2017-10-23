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

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)

#ifndef G2O_IMU
#define G2O_IMU

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "types_sba.h"
#include <Eigen/Geometry>

#define _g 1.1
#define pi 3.14159926
#ifndef BASIC_MATRIX
#define BASIC_MATRIX
Eigen::Matrix3d I3x3 = Eigen::Matrix<double, 3, 3>::Identity();
Eigen::Matrix3d Z3x3 = Eigen::Matrix<double, 3, 3>::Zero();
Eigen::Vector3d I3(1.0, 1.0, 1.0);
Eigen::Vector3d Z3(0.0, 0.0, 0.0);
#endif
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 10, 1> Vector10d;


namespace g2o {
namespace types_IMU {
void init();
}

using namespace Eigen;

class VertexCamera : public BaseVertex<9, Vector9d>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexCamera();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate << Z3, Z3, Z3;
  }

  virtual void oplusImpl(const double* update_)
  {
    SE3Quat temp, temp1;
    temp.fromMinimalVector(estimate().head(6));
    Eigen::Map<const Vector9d> update(update_);

    Vector9d est;
    temp1 = SE3Quat::exp(update.head(6))*temp;      // oprator = ??
    est.head(3)= temp1.rotation().toRotationMatrix().eulerAngles(0, 1, 2);
    est.segment(3,3) = temp1.translation();
    est.tail(3) = estimate().tail(3) + update.tail(3);

    setEstimate(est);
  }
};

class  EdgeIMUpreintegration: public  BaseBinaryEdge<9, Vector9d, VertexCamera, VertexCamera>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::Vector3d g;

  EdgeIMUpreintegration();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  Eigen::Matrix3d Dexp( Eigen::Vector3d theta );
  Eigen::Matrix3d Dlog( Eigen::Vector3d theta );

  void IMUResdual( Vector9d state_i, Vector9d state_j);

  void computeError()   // before compute error, we have to input data from Lib_imupreintegration
  {
    const VertexCamera* v1 = static_cast<const VertexCamera*>(_vertices[0]);
    const VertexCamera* v2 = static_cast<const VertexCamera*>(_vertices[1]);

    IMUResdual( v1->estimate(),  v2->estimate() );
    _error = ResiError;
  }

  virtual void linearizeOplus();

private:
  double delta_tij;
  Eigen::Vector3d delta_pij, delta_vij, bias_a, bias_g, bias_a_previous, bias_g_previous;
  Eigen::Matrix3d delta_Rij;
  Eigen::MatrixXd IMU_cov_ij;     // SIGMA of pre-mea noise

  struct Derivative{
     Eigen::Matrix3d DRij_Dbg;
     Eigen::Matrix3d Dpij_Dba;
     Eigen::Matrix3d Dpij_Dbg;
     Eigen::Matrix3d Dvij_Dba;
     Eigen::Matrix3d Dvij_Dbg;
  }df_dx;

  Eigen::Matrix<double, 9, 1> ResiError;

public:

//  Eigen::Matrix3d Rij_read() { return delta_Rij; }
//  Eigen::Vector3d pij_read() { return delta_pij; }
//  Eigen::Vector3d vij_read() { return delta_vij; }
//  double tij_read() { return delta_tij; }

  Derivative df_dx_read() { return df_dx; }

  void Data_input(Eigen::Matrix3d deltaRij, Eigen::Vector3d deltapij, Eigen::Vector3d deltavij, Eigen::Vector3d biasa, Eigen::Vector3d biasg, Eigen::Vector3d biasa_previous, Eigen::Vector3d biasg_previous, double deltatij);
  void df_dx_input(Eigen::Matrix3d DRij_Dbg, Eigen::Matrix3d Dpij_Dba, Eigen::Matrix3d Dpij_Dbg, Eigen::Matrix3d Dvij_Dba, Eigen::Matrix3d Dvij_Dbg);

};

} // end namespace

#endif
