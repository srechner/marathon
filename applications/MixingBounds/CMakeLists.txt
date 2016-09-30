project(MixingBounds)

set(SOURCE_FILES
        MixingBounds.cpp)

# find required packages
set(BLA_VENDOR Generic)
find_package(ARPACK)
find_package(SuperLU)
find_package(BLAS)

# if all packages could be found
if(ARPACK_FOUND AND SUPERLU_FOUND)

    # compile binary
    add_executable(${PROJECT_NAME} ${SOURCE_FILES})

    # link against libraries
    target_link_libraries(${PROJECT_NAME}
            marathon
            ${ARPACK_LIBRARIES}
            ${SUPERLU_LIBRARIES}
            ${BLAS_LIBRARIES}
            ${MANDATORY_LIBRARIES}
            )

    # install target
    install (TARGETS ${PROJECT_NAME} DESTINATION bin)

else(ARPACK_FOUND AND SUPERLU_FOUND)
    message("-- Warning: Can not build subproject '${PROJECT_NAME}' due to missing dependencies.")
endif(ARPACK_FOUND AND SUPERLU_FOUND)