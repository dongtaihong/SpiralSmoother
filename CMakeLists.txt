cmake_minimum_required(VERSION 2.8)
project(SpiralSmoother)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

# matplotlib-cpp library
find_package(Python3 COMPONENTS NumPy)

# ipopt library
set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/include/coin)
include_directories(
    ${IPOPT_INCLUDE_DIRS} 
    ${CMAKE_SOURCE_DIR}/include
)
link_directories(${IPOPT_LIBRARY_DIRS})
add_library(ipopt_interface ./src/ipopt_interface.cpp)
target_link_libraries(ipopt_interface
                    Python3::Python)               

add_executable(${CMAKE_PROJECT_NAME} 
                src/main.cpp
                src/SpiralSmoother.cpp
)

target_link_libraries(${CMAKE_PROJECT_NAME}  
                        ipopt_interface ipopt
                        Python3::Python
)
