
cmake_minimum_required(VERSION 2.8)

# cmake modules to find 3rd party libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH} )

# Options
option(CUDA "Use CUDA (if available)" ON)

# find 3rd party libraries
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(BLAS REQUIRED)
find_package(ARPACK REQUIRED)
find_package(SuperLU REQUIRED)

# bugfix for Ubuntu 14.04
set(BLAS_LIBRARIES ${BLAS_LIBRARIES} /usr/lib/libblas.so)

# collect common libraries into a single variable
set(MANDATORY_LIBRARIES
        marathon
        ${CMAKE_THREAD_LIBS_INIT}
        ${ARPACK_LIBRARIES}
        ${SUPERLU_LIBRARIES}
        ${BLAS_LIBRARIES}
        )

# find CUDA (optional)
if (CUDA)
    # try to find CUDA
    find_package(CUDA QUIET)
    if (NOT CUDA_FOUND)
        set(CUDA OFF)
    else (NOT CUDA_FOUND)
        set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    endif (NOT CUDA_FOUND)
endif (CUDA)

# Compile Flags
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++11 -O3")

# set include directory for applications (-I)
include_directories("${CMAKE_SOURCE_DIR}/include")

# build library
add_subdirectory(src)

# build applications
add_subdirectory(applications/MixingBounds)
add_subdirectory(applications/RandomSamples)
add_subdirectory(applications/TransitionMatrix)
add_subdirectory(applications/ListStates)
