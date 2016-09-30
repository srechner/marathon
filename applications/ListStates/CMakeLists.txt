project(ListStates)

# compile binary
add_executable(${PROJECT_NAME} ListStates.cpp)


# link against libraries
target_link_libraries(${PROJECT_NAME}
        marathon
    )

# install target
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
