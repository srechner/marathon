cmake_minimum_required(VERSION 3.0)

# Paramters that can be specified to enable optional functions.
# USE_BOOST_SERIALIZATION

# Project Specifics
project(marathon)

# Set Default Build Type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)


# cmake modules to find 3rd party libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# Header Files
set(HEADER_FILES
        ${PROJECT_SOURCE_DIR}/include/marathon/spectral_bound.h
        ${PROJECT_SOURCE_DIR}/include/marathon/markov_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/path_congestion.h
        ${PROJECT_SOURCE_DIR}/include/marathon/path_construction_scheme.h
        ${PROJECT_SOURCE_DIR}/include/marathon/basic_random.h
        ${PROJECT_SOURCE_DIR}/include/marathon/rational.h
        ${PROJECT_SOURCE_DIR}/include/marathon/state.h
        ${PROJECT_SOURCE_DIR}/include/marathon/count.h
        ${PROJECT_SOURCE_DIR}/include/marathon/enumerate.h
        ${PROJECT_SOURCE_DIR}/include/marathon/statistics.h
        ${PROJECT_SOURCE_DIR}/include/marathon/state_graph.h
        ${PROJECT_SOURCE_DIR}/include/marathon/transition.h
        ${PROJECT_SOURCE_DIR}/include/marathon/transition_matrix.h
        ${PROJECT_SOURCE_DIR}/include/marathon/integer.h
        ${PROJECT_SOURCE_DIR}/include/marathon/mixing_time.h
        ${PROJECT_SOURCE_DIR}/include/marathon/structure.h
        ${PROJECT_SOURCE_DIR}/include/marathon/combination_generator.h
        ${PROJECT_SOURCE_DIR}/include/marathon/permutation_generator.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/common.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/binary_matrix.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/random_generator.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/metrics.h
        ${PROJECT_SOURCE_DIR}/include/marathon/random_generator.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/count.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/decompose.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/markov_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/random_generator_exact.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/random_generator_intelligent.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/random_generator_mcmc.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/instance.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/curveball.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/canonical_path_kannan.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/switch_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/switch_chain_canonical.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/edge_switch_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/fixed_margin/switch_chain_intelligent.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/count.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/enumerate.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/markov_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/random_generator_mcmc.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/random_generator_exact.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/realize.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/instance.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/informed_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/informed_chain_dynamic.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/simple_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/binary_matrix/interval_margin/simple_chain_dynamic.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/bipartite_matching.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/markov_chain.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/count.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/enumerate.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/realize.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/metrics.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/random_generator_mcmc.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/random_generator_exact.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/perfect_matching_broder86.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/canonical_path_js89.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/perfect_matching_jsv04.h
        ${PROJECT_SOURCE_DIR}/include/marathon/matching/sparse_bipartite_graph.h
        )

# add marathon include path
include_directories(${PROJECT_SOURCE_DIR}/include)

# Compile Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# find boost headers
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# find boost serialization
if(USE_BOOST_SERIALIZATION)
    find_package(Boost REQUIRED COMPONENTS serialization)
    add_definitions(-DUSE_BOOST_SERIALIZATION)
    include_directories(${Boost_INCLUDE_DIRS})
    set(MARATHON_LIBRARIES ${MARATHON_LIBRARIES} ${Boost_LIBRARIES})
endif()

# try to find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

# try to find blas for linear algebra
set(BLA_VENDOR Generic)
find_package(BLAS)
if(BLAS_FOUND)
    add_definitions(-DUSE_BLAS)
    message(STATUS "Found BLAS libraries ${BLAS_LIBRARIES}")
    include_directories(${BLAS_INCLUDE_DIRS})
    set(MARATHON_LIBRARIES ${MARATHON_LIBRARIES} ${BLAS_LIBRARIES})
endif()

# try to find armadillo for linear algebra
find_package(Armadillo)
if(ARMADILLO_FOUND)
    add_definitions(-DUSE_ARMADILLO)
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    set(MARATHON_LIBRARIES ${MARATHON_LIBRARIES} ${ARMADILLO_LIBRARIES})
endif()

# try to find Eigen3
find_package(Eigen3)
if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
    set(MARATHON_LIBRARIES ${MARATHON_LIBRARIES} ${EIGEN3_LIBRARIES})
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Using Libraries: ${MARATHON_LIBRARIES}")

# Install Target
install(FILES ${HEADER_FILES} DESTINATION include/marathon/)

# build examples
add_subdirectory(examples/sample)
add_subdirectory(examples/count)
add_subdirectory(examples/enumerate)
add_subdirectory(examples/transitionMatrix)
add_subdirectory(examples/totalMixingTime)
add_subdirectory(examples/spectralBound)
