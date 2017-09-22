cmake_minimum_required(VERSION 2.8)

project(Imupreintegration_Test)
set( CMAKE_BUILD_TYPE "Debug" )

INCLUDE_DIRECTORIES( "/usr/include/eigen3" )    # Eigen

# 添加自己写的imu preintegration
add_library( preintegration "Lib_preintegration.cpp" )
add_library( preintegration_shared SHARED Lib_preintegration.cpp )

# 主函数可执行文件
add_executable(${PROJECT_NAME} "main.cpp")

# 将库文件链接到可执行程序上
target_link_libraries( ${PROJECT_NAME} preintegration_shared )
