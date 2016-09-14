project(RandomSamples)

# compile binary
add_executable(${PROJECT_NAME} RandomSamples.cpp)

# link against libraries
target_link_libraries(${PROJECT_NAME} ${MANDATORY_LIBRARIES})

# install target
install (TARGETS ${PROJECT_NAME} DESTINATION bin)