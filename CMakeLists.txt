cmake_minimum_required(VERSION 2.8)

project(Imupreintegration_Test)
set( CMAKE_BUILD_TYPE "Debug" )

#find_package( G2O REQUIRED )

find_package(Eigen3)
include_directories(
${PROJECT_SOURCE_DIR}
${EIGEN3_INCLUDE_DIR}
)
# 添加自己写的imu preintegration
add_library( preintegration "Lib_preintegration.cpp" )
add_library( preintegration_shared SHARED Lib_preintegration.cpp )

# 添加优化
add_library( optimizeIMU "optimizeIMU.cpp" )
add_library( optimizeIMU_shared SHARED optimizeIMU.cpp )

# 主函数可执行文件
add_executable(${PROJECT_NAME} "main.cpp")

# 将库文件链接到可执行程序上
target_link_libraries( ${PROJECT_NAME} 
preintegration_shared
optimizeIMU_shared
${PROJECT_SOURCE_DIR}/g2o/lib/libg2o.so
)
