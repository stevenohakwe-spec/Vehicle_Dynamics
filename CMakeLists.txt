
cmake_minimum_required(VERSION 3.21)

project(VehicleDynamicsSimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SFML 3 REQUIRED COMPONENTS Graphics)

add_executable(VehicleDynamics main.cpp)

target_link_libraries(
    VehicleDynamics
    PRIVATE
    SFML::Graphics
)
