# Project Specifics
project(marathon)

# version number
set(marathon_VERSION_MAJOR 0)
set(marathon_VERSION_MINOR 4)

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_SOURCE_DIR}/include/marathon/config.h.in"
        "${CMAKE_SOURCE_DIR}/include/marathon/config.h"
)

# Include Files
set(HEADER_FILES
        ../include/marathon/Combinatorics.h
        ../include/marathon/config.h
        ../include/marathon/CudaWrapper.h
        ../include/marathon/Decimal.h
        ../include/marathon/Eigenvalues.h
        ../include/marathon/marathon.h
        ../include/marathon/MarkovChain.h
        ../include/marathon/PathCongestion.h
        ../include/marathon/PathConstructionScheme.h
        ../include/marathon/Random.h
        ../include/marathon/Rational.h
        ../include/marathon/State.h
        ../include/marathon/StateGraph.h
        ../include/marathon/Transition.h
        ../include/marathon/TransitionMatrix.h
        ../include/marathon/TransitionMatrixCBLAS.h
        ../include/marathon/TransitionMatrixCuBLAS.h
        ../include/marathon/TransitionMatrixCuBLASXt.h
        ../include/marathon/chain/Chains.h
        ../include/marathon/chain/bipgraph/BinaryMatrix.h
        ../include/marathon/chain/bipgraph/Curveball.h
        ../include/marathon/chain/bipgraph/HavelHakimi.h
        ../include/marathon/chain/bipgraph/KannanCanPath.h
        ../include/marathon/chain/bipgraph/SwitchChain.h
        ../include/marathon/chain/bipgraph/SwitchChainBerger.h
        ../include/marathon/chain/matching/BipartiteMatching.h
        ../include/marathon/chain/matching/Broder86.h
        ../include/marathon/chain/matching/JS89CanPath.h
        ../include/marathon/chain/matching/JSV04.h
        ../include/marathon/chain/matching/SparseBipartiteGraph.h
        )

# C++ Source Files
set(CPP_SOURCE_FILES
        marathon/marathon.cpp
        marathon/Rational.cpp
        marathon/StateGraph.cpp
        marathon/Transition.cpp
        marathon/MarkovChain.cpp
        marathon/PathCongestion.cpp
        marathon/Eigenvalues.cpp
        marathon/TotalMixingTime.cpp
        marathon/ShortestPaths.cpp
        marathon/InitFinalize.cpp
        marathon/Random.cpp
        marathon/Combinatorics.cpp
        marathon/TransitionMatrixCBLAS.cpp
        marathon/TransitionMatrixCuBLAS.cpp
        marathon/TransitionMatrixCuBLASXt.cpp
        marathon/chain/matching/BipartiteMatching.cpp
        marathon/chain/matching/Broder86.cpp
        marathon/chain/matching/JSV04.cpp
        marathon/chain/matching/SparseBipartiteGraph.cpp
        marathon/chain/matching/JS89CanPath.cpp
        marathon/chain/bipgraph/Curveball.cpp
        marathon/chain/bipgraph/SwitchChain.cpp
        marathon/chain/bipgraph/SwitchChainBerger.cpp
        marathon/chain/bipgraph/BinaryMatrix.cpp
        marathon/chain/bipgraph/HavelHakimi.cpp
        marathon/chain/bipgraph/KannanCanPath.cpp
        )

# CUDA Source Files
set(CUDA_SOURCE_FILES
        marathon/CudaVariationDistance.cu
        marathon/CudaWrapper.cu
        marathon/CudaInitFinalize.cu
        )

# compile CUDA files
if (CUDA)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)                      # turn off flag propagation
    set(EXTRA_INCLUDES ${CUDA_INCLUDE_DIRS})                # use cuda headers
    set(EXTRA_LIBS ${EXTRA_LIBS} ${CUDA_LIBRARIES})         # use cuda libraries
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}                  # set nvcc compile flags
            -std=c++11 -O3
            -Xcompiler -fopenmp
            -Xcompiler -fPIC
            -D_FORCE_INLINES
            -D_MWAITXINTRIN_H_INCLUDED)
    cuda_compile(cuda_objs ${CUDA_SOURCE_FILES})            # compile .cu sources
endif (CUDA)

# Compile library
add_library(marathon STATIC ${HEADER_FILES} ${CPP_SOURCE_FILES} ${cuda_objs})

# Install
install (TARGETS marathon DESTINATION lib)
install (FILES ${HEADER_FILES} DESTINATION include/marathon/)
