macro(ADD_EXAMPLE name)
    add_executable(example_${name} example_${name}.f90)
    target_link_libraries(example_${name} "${PROJECT_NAME}")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:example_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADD_EXAMPLE)

add_subdirectory(mpc)
add_subdirectory(compact_mpc)
add_subdirectory(markowitz)

