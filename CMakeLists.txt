cmake_minimum_required(VERSION 3.16)
project(Num_methods)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS " -fconcepts -lboost_unit_test_framework-gcc-1_61")# -Weffc++ -pedantic -Wextra -Wall -Wnon-virtual-dtor -Wshadow -Winit-self -Wredundant-decls -Wcast-align -Wundef -Wfloat-equal -Winline -Wunreachable-code -Wmissing-include-dirs -Wswitch-enum -Wall -Wctor-dtor-privacy -Wold-style-cast -Woverloaded-virtual -Wsign-promo -Wduplicated-branches -Wduplicated-cond -Wcast-qual -Wconversion -Wsign-conversion -Wlogical-op")

include_directories(Lib)

add_executable(Num_methods
        Lib/Math.cpp
        Lib/Matrix.cpp
        Run_tests.cpp
        )

include(PVS-Studio.cmake)
pvs_studio_add_target(TARGET pvs-analyze ALL
                      OUTPUT FORMAT errorfile
                      ANALYZE ${PROJECT_NAME}
                      MODE GA:1,2 OP
                      LOG target.err)