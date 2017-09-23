#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Lib_preintegration.h"
#include "se3.h"

using namespace std;

int main(int argc, char *argv[])
{

    ImuPreintegration IMUP;

/// ROS bag data reading
///
    ifstream Imudata;
    int flag_acc = 15, flag_gyro = 27, data_cont = 0;
    vector<double> simul_meas, time_seq;    // save data during a measurement
    vector< vector<double> > Imu_meas;
    string IMUFilePath = "/home/gentle/Documents/Imu_preintegration/Imupreintegration_Test/IMUDataRawData.txt";
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
                else if ( (data_cont >= flag_acc && data_cont <= flag_acc+2) || (data_cont >= flag_gyro && data_cont <= flag_gyro+2))
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
    Eigen::Vector3d acc, gyro;
    for ( vector< vector<double> >::iterator it = Imu_meas.begin(); it != Imu_meas.end(); ++it  )
    {
        simul_meas = *it;
        data_cont = 0;
        for ( vector<double>::iterator it_elem = simul_meas.begin(); it_elem != simul_meas.end(); ++it_elem )
        {
            if ( t_cont < time_seq.size()-2 ){
                dt = ( time_seq[t_cont+1] - time_seq[t_cont] ) / 1e9 ;    // 最后一个时刻沿用前一时刻dt
            }

            if(data_cont < 3)
            {
                gyro(data_cont) = *it_elem;
            }
            else
            {
                acc(data_cont-3) = *it_elem;
            }

            data_cont += 1;
            t_cont += 1;
        }

        IMUP.Preintegration(acc, gyro, dt);

    }

    return 0;
}
