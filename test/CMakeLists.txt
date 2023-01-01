set(ProjectTest ${CMAKE_PROJECT_NAME}_Tests)

set(Sources
    ExampleTests.cpp
)

add_executable(${ProjectTest} ${Sources})

target_link_libraries(
    ${ProjectTest}
    PUBLIC ${CMAKE_PROJECT_NAME}_lib
    gtest_main
)

add_test(
    NAME ${ProjectTest}
    COMMAND ${ProjectTest}
)

