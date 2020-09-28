cmake_minimum_required(VERSION 3.17)
project(GBQ)

set(CMAKE_CXX_STANDARD 14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# CPLEX
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD")

# GUROBI
#set(GRB_LICENSE_FILE /home/carlosrey/gurobi903/linux64/gurobi.lic)
set(GUROBI_DIR /home/carlosrey/gurobi903/linux64)
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

# OPEN MP
find_package(OpenMP REQUIRED)

add_executable(GBQ main.cpp localsearch.cpp localsearch.h OPTUtils.h utils.h utils.cpp Solution.cpp Solution.h QMKP_01qpCPX2.cpp QMKP_01qpCPX2.h qkp_grdy.cpp qkp_grdy.h quadknap.c)

message("${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} ${GUROBI_CXX_DEBUG_LIBRARY}")
target_link_libraries(GBQ ${CPLEX_LIBRARIES} pthread)

target_link_libraries(GBQ optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
target_link_libraries(GBQ  ${GUROBI_LIBRARY})
target_link_libraries(GBQ OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()