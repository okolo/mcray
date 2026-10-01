add_executable(table_function_log_zeros table_function_log_zeros.cpp)
target_include_directories(table_function_log_zeros PRIVATE ../src/lib)
target_link_libraries(table_function_log_zeros mcray)
add_test(NAME table_function_log_zeros COMMAND table_function_log_zeros)

add_executable(psb_branching_bounds psb_branching_bounds.cpp)
target_include_directories(psb_branching_bounds PRIVATE ../src/lib)
target_link_libraries(psb_branching_bounds mcray)
add_test(NAME psb_branching_bounds COMMAND psb_branching_bounds)

add_test(NAME crbeam_help COMMAND CRbeam --help)
add_test(NAME z2t_smoke COMMAND z2t 0.1)
set_tests_properties(table_function_log_zeros psb_branching_bounds crbeam_help z2t_smoke
    PROPERTIES TIMEOUT 30 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin")
