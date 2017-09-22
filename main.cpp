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
//    ImuPreintegration IMUP;

/// ROS bag data reading
    ifstream Imudata;
    vector<double> single_line;
    vector< vector<double> > Imu_meas;
    string IMUFilePath = "/home/gentle/Documents/Imu_preintegration/Imupreintegration_Test/IMUDataData.txt";
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
        if(!s.empty())
        {
            ss << s;
            while ( ss ) {
                ss >> data;
                single_line.push_back(data);
            }

        }
        Imu_meas.push_back(single_line);
        single_line.clear();
    }

/* preintegration
 */
    //IMUP.Preintegration(acc, gyro, dt);

    return 0;
}
