#include <g2o/g2o/core/block_solver.h>
#include <g2o/g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/g2o/solvers/linear_solver_eigen.h>
#include <g2o/g2o/core/robust_kernel_impl.h>
#include <g2o/g2o/solvers/linear_solver_dense.h>
#include <g2o/g2o/types/types_IMU.h>

#include "Lib_preintegration.h"
#include "se3.h"

void optimizeIMUPre( ImuPreintegration& origin_IMUP,  Obj_State& state_i, Obj_State& state_j)
{
    // 构造优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_9_9::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_9_9::PoseMatrixType>();

    g2o::BlockSolver_9_9 * solver_ptr = new g2o::BlockSolver_9_9(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // 添加位姿结点
    g2o::VertexCamera* vSE3 = new g2o::VertexCamera();
    vSE3->setEstimate(state_i.state_j);
    vSE3->setId(state_i.Id);
    vSE3->setFixed(state_i.Id == 0);   // 第一个结点是否fix
    optimizer.addVertex(vSE3);

    // 添加位姿结点
    g2o::VertexCamera* vSE3_j = new g2o::VertexCamera();
    vSE3_j->setEstimate(state_j.state_j);
    vSE3_j->setId(state_j.Id);
    vSE3_j->setFixed(false);      // ## 第一个结点是否fix住？
    optimizer.addVertex(vSE3_j);


    //添加IMUEdge
    SE3Quat IMUP_temp(origin_IMUP.Rij_read(), origin_IMUP.pij_read());
    Eigen::Matrix<double,9,1> obs;
    obs << IMUP_temp.log(), origin_IMUP.vij_read();

    g2o::EdgeIMUpreintegration* e = new g2o::EdgeIMUpreintegration();

    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(state_i.Id)));
    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(state_j.Id)));
    e->setMeasurement(obs);

    e->setInformation(origin_IMUP.IMUCov_read());

    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
    e->setRobustKernel(rk);
    rk->setDelta(sqrt(10));     // ## 10 is not reliable

    e->Data_input( origin_IMUP.Rij_read(), origin_IMUP.pij_read(), origin_IMUP.vij_read(), origin_IMUP.biasa_read(), origin_IMUP.biasg_read(), origin_IMUP.bias_prea_read(), origin_IMUP.bias_preg_read(), origin_IMUP.tij_read());
    e->df_dx_input(origin_IMUP.df_dx_read().DRij_Dbg, origin_IMUP.df_dx_read().Dpij_Dba, origin_IMUP.df_dx_read().Dpij_Dbg, origin_IMUP.df_dx_read().Dvij_Dba, origin_IMUP.df_dx_read().Dvij_Dbg);

    optimizer.addEdge(e);

//    vpEdgesStereo.push_back(e);
//    vnIndexEdgeStereo.push_back(i);

    const int its[4]={10,10,10,10};     // 四次迭代，每次迭代的次数

    for(size_t it=0; it<4; it++)
    {
//        vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));    // 之前有加锁，所以再赋值
        optimizer.initializeOptimization(0);    // 对level为0的边进行优化
        optimizer.optimize(its[it]);

//        for(size_t i=0, iend = vpEdgesStereo.size(); i<iend; i++)
//        {
//            g2o::EdgeIMUpreintegration* e = vpEdgesStereo[i];

            e->computeError();
            e->setLevel(0);

            if(it==2)
                e->setRobustKernel(0);
//        }

        if(optimizer.edges().size()<10)
            break;
    }

    g2o::VertexCamera* camera_recov = static_cast<g2o::VertexCamera*>(optimizer.vertex(state_i.Id));
    Vector9d statei_recov = camera_recov->estimate();
    state_i = statei_recov;

    camera_recov = static_cast<g2o::VertexCamera*>(optimizer.vertex(state_j.Id));
    Vector9d statej_recov = camera_recov->estimate();
    state_j = statej_recov;

}
