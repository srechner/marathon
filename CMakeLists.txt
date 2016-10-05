cmake_minimum_required(VERSION 2.8)

# cmake modules to find 3rd party libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# Project Specifics
project(marathon)

# version number
set(marathon_VERSION_MAJOR 0)
set(marathon_VERSION_MINOR 4)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/include/marathon/config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
)

# Include Files
set(HEADER_FILES
        include/marathon/Combinatorics.h
        include/marathon/Eigenvalue.h
        include/marathon/marathon.h
        include/marathon/MarkovChain.h
        include/marathon/PathCongestion.h
        include/marathon/PathConstructionScheme.h
        include/marathon/Random.h
        include/marathon/Rational.h
        include/marathon/State.h
        include/marathon/StateGraph.h
        include/marathon/Transition.h
        include/marathon/TransitionMatrix.h
        include/marathon/TransitionMatrixCBLAS.h
        include/marathon/TransitionMatrixCuBLAS.h
        include/marathon/TransitionMatrixCuBLASXt.h
        include/marathon/bipgraph/BinaryMatrix.h
        include/marathon/bipgraph/HavelHakimi.h
        include/marathon/bipgraph/KannanCanPath.h
        include/marathon/bipgraph/SwitchChain.h
        include/marathon/bipgraph/SwitchChainBerger.h
        include/marathon/matching/BipartiteMatching.h
        include/marathon/matching/Broder86.h
        include/marathon/matching/JS89CanPath.h
        include/marathon/matching/JSV04.h
        include/marathon/matching/SparseBipartiteGraph.h
        include/marathon/Types.h
        include/marathon/MixingTime.h
        include/marathon/Diameter.h
        include/marathon/Memory.h
        )

# CUDA Source Files
set(CU_SOURCE_FILES
        src/marathon/cuda/CudaVariationDistance.cu)

# find mandatory libraries
set(BLA_VENDOR Generic)
find_package(Boost REQUIRED)
find_package(Threads)
find_package(OpenMP)
find_package(CUDA)
find_package(ARPACK)
find_package(SuperLU)
find_package(BLAS)
find_package(Eigen3)

# Compile Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# set include directory for applications (-I)
include_directories(
        ${PROJECT_SOURCE_DIR}/include   # marathon/include
        ${PROJECT_BINARY_DIR}           # local build folder (where config.h is located)
        ${Boost_INCLUDE_DIRS}
)

# if OpenMP was found
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

#if a thread library was found
if(Threads_FOUND)
    set(BASIC_LIBRARIES ${BASIC_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})     # add thread library
endif(Threads_FOUND)

# if CUDA is located at the system
if (CUDA_FOUND)
    add_definitions(-DCUDA)                                                     # enable CUDA flag
    cuda_add_library(${PROJECT_NAME} ${CU_SOURCE_FILES})        # compile .cu sources
    set(BASIC_LIBRARIES ${BASIC_LIBRARIES}
            marathon ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})                # add cuda libraries
    install(TARGETS marathon DESTINATION lib)                                   # Install Target
endif (CUDA_FOUND)

# Install Target
install(FILES ${HEADER_FILES} ${PROJECT_BINARY_DIR}/config.h DESTINATION include/marathon/)

# build applications
add_subdirectory(src/apps/MixingBounds)
add_subdirectory(src/apps/TransitionMatrix)
add_subdirectory(src/apps/ListStates)
add_subdirectory(src/apps/RandomSamples)
