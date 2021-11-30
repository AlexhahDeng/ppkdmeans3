cmake_minimum_required (VERSION 3.13)
set(CMAKE_CXX_STANDARD 17)

project (dyx)

include_directories("./install")

set(SCI_DIR /home/dd/Ubu20/test/install/lib/cmake/SCI)
find_package(SCI REQUIRED)


add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} SCI::SCI-OT SCI::SCI-HE)
