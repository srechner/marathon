project(TransitionMatrix)

set(SOURCE_FILES
        TransitionMatrix.cpp)

# if all packages could be found
if (BLAS_FOUND AND Boost_FOUND)

    # compile binary
    add_executable(${PROJECT_NAME} ${SOURCE_FILES})

    # link against libraries
    target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARIES} ${BASIC_LIBRARIES})

    # install target
    install(TARGETS ${PROJECT_NAME} DESTINATION bin)

else (BLAS_FOUND AND Boost_FOUND)

    message("-- Warning: Can not build subproject '${PROJECT_NAME}' due to missing dependencies.")

endif (BLAS_FOUND AND Boost_FOUND)