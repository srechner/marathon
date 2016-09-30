project(RandomSamples)

set(SOURCE_FILES
        RandomSamples.cpp)

# find required packages

# if all packages could be found
if(TRUE)

    # compile binary
    add_executable(${PROJECT_NAME} ${SOURCE_FILES})

    # link against libraries
    target_link_libraries(${PROJECT_NAME}
            marathon
            ${MANDATORY_LIBRARIES}
            )

    # install target
    install (TARGETS ${PROJECT_NAME} DESTINATION bin)

else(TRUE)
    message("-- Warning: Can not build subproject '${PROJECT_NAME}' due to missing dependencies.")
endif(TRUE)