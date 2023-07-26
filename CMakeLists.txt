cmake_minimum_required (VERSION 3.22)
project(A)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the path to where the RDKit installed
set(RDKit_DIR "/home/luigi/mambaforge/envs/my-rdkit-env/lib/cmake/rdkit")
find_package(RDKit REQUIRED)

set(BOOST_ROOT "/home/luigi/boost")
find_package(Boost COMPONENTS timer system REQUIRED)

# Define the executable
add_executable(main main.cpp)

target_link_libraries(main RDKit::SmilesParse RDKit::MolHash Boost::timer)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BOOST_ROOT "/ home/luigi/boost")


find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS timer system REQUIRED)
target_link_libraries(main RDKit::GraphMol  Boost::timer)
