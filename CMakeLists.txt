cmake_minimum_required(VERSION 3.16)
project(Num_methods)

set(CMAKE_CXX_STANDARD 14)

include_directories(2_Lab/Lib1)
include_directories(3_Lab)
include_directories(Lib)

add_executable(Num_methods
        1_Lab/main.cpp
        2_Lab/Lib1/Different.cpp
        2_Lab/Lib1/Different.h
        2_Lab/Lib1/List.cpp
        2_Lab/Lib1/List.h
        2_Lab/Lib1/Matrix.cpp
        2_Lab/Lib1/Matrix.h
        2_Lab/Lib1/Tree.h
        2_Lab/Lib1/Vector.cpp
        2_Lab/Lib1/Vector.h
        2_Lab/main.cpp
        3_Lab/main.cpp
        3_Lab/version.h
        4_Lab/main.cpp
        5_Lab/main.cpp
        6_Lab/main.cpp
        6_Lab_default/main.cpp
        7_Lab/main.cpp
        7_Lab/main_default.cpp
        7_Lab_default/main.cpp
        Lab_6-7/7_Lab_default/main.cpp
        Lib/Different.cpp
        Lib/Different.h
        Lib/List.cpp
        Lib/List.h
        Lib/Math.cpp
        Lib/Math.h
        Lib/Matrix.cpp
        Lib/Matrix.h
        Lib/Tree.h
        Lib/Vector.cpp
        Lib/Vector.h)
