cmake_minimum_required(VERSION 3.9)

add_executable(main main.cpp individual.cpp genetic_algorithm.cpp random.cpp)
set_property(TARGET main PROPERTY CXX_STANDARD 17)