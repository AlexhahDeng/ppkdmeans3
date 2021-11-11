cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(pp3
       	LANGUAGES CXX)

find_package(helib REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})

# add_executable(comparison_circuit comparison_circuit.cpp comparator.cpp tools.cpp)

add_executable(pp3 my_tools.cpp cloud.cpp comparator.cpp func.cpp pp3.cpp tools.cpp)
# add_executable(test_sc cloud.cpp comparator.cpp func.cpp test_sc.cpp tools.cpp)


target_link_libraries(pp3 helib)


# cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

# ## Use -std=c++17 as default.
# set(CMAKE_CXX_STANDARD 17)
# ## Disable C++ extensions
# set(CMAKE_CXX_EXTENSIONS OFF)
# ## Require full C++ standard
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# project(test_sc
#        	LANGUAGES CXX)

# find_package(helib REQUIRED)

# include_directories(${PROJECT_SOURCE_DIR})

# # add_executable(comparison_circuit comparison_circuit.cpp comparator.cpp tools.cpp)

# # add_executable(pp3 my_tools.cpp cloud.cpp comparator.cpp func.cpp pp3.cpp tools.cpp)
# add_executable(test_sc my_tools.cpp cloud.cpp comparator.cpp func.cpp test_sc.cpp tools.cpp)

# target_link_libraries(test_sc helib)