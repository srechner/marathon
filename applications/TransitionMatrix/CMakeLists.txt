project(TransitionMatrix)

set(SOURCE_FILES
        TransitionMatrix.cpp)

# find required packages
set(BLA_VENDOR Generic)
find_package(BLAS)

# if all packages could be found
if(BLAS_FOUND)

    # compile binary
    add_executable(${PROJECT_NAME} ${SOURCE_FILES})

    # link against libraries
    target_link_libraries(${PROJECT_NAME}
            marathon
            ${BLAS_LIBRARIES}
            ${MANDATORY_LIBRARIES}
            )

    # install target
    install (TARGETS ${PROJECT_NAME} DESTINATION bin)

else(BLAS_FOUND)
    message("-- Warning: Can not build subproject '${PROJECT_NAME}' due to missing dependencies.")
endif(BLAS_FOUND)