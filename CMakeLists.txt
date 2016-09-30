cmake_minimum_required(VERSION 2.8)

# cmake modules to find 3rd party libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# Project Specifics
project(marathon)

# version number
set(marathon_VERSION_MAJOR 0)
set(marathon_VERSION_MINOR 4)

# Include Files
set(HEADER_FILES
        include/marathon/Combinatorics.h
        include/marathon/Cuda.h
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
        include/marathon/chain/Chains.h
        include/marathon/chain/bipgraph/BinaryMatrix.h
        include/marathon/chain/bipgraph/HavelHakimi.h
        include/marathon/chain/bipgraph/KannanCanPath.h
        include/marathon/chain/bipgraph/SwitchChain.h
        include/marathon/chain/bipgraph/SwitchChainBerger.h
        include/marathon/chain/matching/BipartiteMatching.h
        include/marathon/chain/matching/Broder86.h
        include/marathon/chain/matching/JS89CanPath.h
        include/marathon/chain/matching/JSV04.h
        include/marathon/chain/matching/SparseBipartiteGraph.h
        include/marathon/Types.h
        include/marathon/MixingTime.h
        include/marathon/Diameter.h
        include/marathon/Memory.h
        )

# C++ Source Files
set(CPP_SOURCE_FILES
        src/marathon/marathon.cpp
        src/marathon/Rational.cpp
        src/marathon/StateGraph.cpp
        src/marathon/Transition.cpp
        src/marathon/MarkovChain.cpp
        src/marathon/PathCongestion.cpp
        src/marathon/Eigenvalue.cpp
        src/marathon/ShortestPaths.cpp
        src/marathon/Random.cpp
        src/marathon/Combinatorics.cpp
        src/marathon/TransitionMatrixCBLAS.cpp
        src/marathon/TransitionMatrixCuBLAS.cpp
        src/marathon/TransitionMatrixCuBLASXt.cpp
        src/marathon/chain/matching/BipartiteMatching.cpp
        src/marathon/chain/matching/Broder86.cpp
        src/marathon/chain/matching/JSV04.cpp
        src/marathon/chain/matching/SparseBipartiteGraph.cpp
        src/marathon/chain/matching/JS89CanPath.cpp
        src/marathon/chain/bipgraph/SwitchChain.cpp
        src/marathon/chain/bipgraph/SwitchChainBerger.cpp
        src/marathon/chain/bipgraph/BinaryMatrix.cpp
        src/marathon/chain/bipgraph/HavelHakimi.cpp
        src/marathon/chain/bipgraph/KannanCanPath.cpp
        src/marathon/Cuda.cpp
        )

# CUDA Source Files
set(CU_SOURCE_FILES
        src/marathon/CudaVariationDistance.cu
        src/marathon/CudaWrapper.cu
        src/marathon/CudaInitFinalize.cu
        )

# check whether to find a CUDA installation
find_package(CUDA)
if (CUDA_FOUND)
    set(CUDA ON)        # enable CUDA flag
    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CUDA_LIBRARIES} ${CUBLAS_LIBRARIES})
endif (CUDA_FOUND)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/include/marathon/config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
)

# set include directory for applications (-I)
include_directories(
        ${PROJECT_SOURCE_DIR}/include   # marathon/include
        ${PROJECT_BINARY_DIR}           # local build folder (where config.h is located)
)

# Compile Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")

# if CUDA is located at the system
if (CUDA_FOUND)
    set(MANDATORY_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    include_directories(${CUDA_INCLUDE_DIRS})
    cuda_compile(CUDA_OBJECTS ${CU_SOURCE_FILES})            # compile .cu sources
endif (CUDA_FOUND)

# Compile library
add_library(marathon STATIC ${HEADER_FILES} ${CPP_SOURCE_FILES} ${CUDA_OBJECTS})

# Install
install (TARGETS marathon DESTINATION lib)
install (FILES ${HEADER_FILES} ${PROJECT_BINARY_DIR}/config.h DESTINATION include/marathon/)

# build applications
add_subdirectory(applications/MixingBounds)
add_subdirectory(applications/RandomSamples)
add_subdirectory(applications/TransitionMatrix)
add_subdirectory(applications/ListStates)