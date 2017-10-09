#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Lib_preintegration.h"
#include "se3.h"

using namespace std;


/// used to record state
struct Obj_State {
    Eigen::Matrix3d R;
    Eigen::Vector3d p;
    Eigen::Vector3d v;

    Vector9d state_j;

};


int main(int argc, char *argv[])
{
    int update_step = 50;
    Obj_State State_i;
    ImuPreintegration IMUP;

    State_i.R = Eigen::Matrix<double, 3, 3>::Identity();
    State_i.p << 0, 0, 0;
    State_i.v << 0, 0, 0;

/// ROS bag data reading
///
    ifstream Imudata;
    int flag_AngVel = 15, flag_LinAcc = 27, data_cont = 0;
    vector<double> simul_meas, time_seq;    // save data during a measurement
    vector< vector<double> > Imu_meas;
    string IMUFilePath = "/home/gentle/Documents/Imu_preintegration/Imupreintegration_Test/IMUcam_RawData.txt";
    Imudata.open( IMUFilePath.c_str() );
    while(!Imudata.eof())
    {
        string s;
        stringstream ss;
        double data = 0;

        //将字符串中的','用空格代替
        getline(Imudata,s);
        int pos = s.find(',');
        while (pos != string::npos)
        {
            s = s.replace(pos, 1, 1, ' ');
            pos = s.find(',');
        }

        // 读入数据
        data_cont = 0;
        if(!s.empty())
        {
            ss << s;
            while ( ss ) {

                ss >> data;
                if (data_cont == 1)     // save time stamps
                {
                    time_seq.push_back(data);
                }
                else if ( (data_cont >= flag_AngVel && data_cont <= flag_AngVel+2) || (data_cont >= flag_LinAcc && data_cont <= flag_LinAcc+2))
                {
                    simul_meas.push_back(data);
                }
                data_cont += 1;
            }

        }
        Imu_meas.push_back(simul_meas);
        simul_meas.clear();
    }


/// Preintegration
///
    int t_cont = 0;
    double dt;
    SE3Quat exp_temp;
    Vector6d r_t;
    Eigen::Vector3d acc, gyro;
    ofstream f_save("ImuPre_result.txt", ios::app);  // 结果存储文件
    for ( vector< vector<double> >::iterator it = Imu_meas.begin(); it != Imu_meas.end(); ++it  )
    {
        simul_meas = *it;
        data_cont = 0;

        // data extraction
        for ( vector<double>::iterator it_elem = simul_meas.begin(); it_elem != simul_meas.end(); ++it_elem )
        {
            if(data_cont < 3)
            {
                gyro(data_cont) = *it_elem;     //  x,y,z rotation velocity, coordinate: ?  (9250而言 朝向旋转轴方向y,z顺时针为负，逆时针方向为正, x相反)
            }
            else
            {
                acc(data_cont-3) = *it_elem;    //  x,y,z  coordinate: ENU or ?
            }

            data_cont += 1;
        }
        if ( t_cont < time_seq.size()-2 )
        {
            dt = ( time_seq[t_cont+1] - time_seq[t_cont] ) / 1e9 ;    // 最后一个时刻沿用前一时刻dt
            t_cont += 1;
        }

        IMUP.Preintegration(acc, gyro, dt);

        if ( t_cont % update_step == 0)
        {
            State_i.state_j = IMUP.predict(State_i.R, State_i.p, State_i.v);
            IMUP.Reset();

            r_t << State_i.state_j.head(3), 0,0,0;
            State_i.R = exp_temp.exp(r_t).rotation().matrix();
            State_i.p = State_i.state_j.segment(3,3);
            State_i.v = State_i.state_j.tail(3);

            if (!f_save)
            {
                cout<<"failed to open file."<<endl;
                exit (0);
            }
            else
            {
                f_save << fixed << setprecision(8) << State_i.state_j.transpose() << endl;
            }
        }
    }

    f_save.close();
    return 0;
}
