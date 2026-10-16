# Also usable by the portable server overlay on the fetched revision.
if(NOT TARGET em_secondary_probe)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    set(sampling_generated "${CMAKE_CURRENT_BINARY_DIR}/sampling_reference")
    execute_process(COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/sampling/prepare_reference.py"
        --repo "${CMAKE_CURRENT_SOURCE_DIR}/.." --nr "${MCRAY_NR_REFERENCE_DIR}"
        --output "${sampling_generated}" COMMAND_ERROR_IS_FATAL ANY)
    add_executable(em_secondary_probe sampling/secondary_probe.cpp
        "${sampling_generated}/historical_sampler.cpp"
        "${sampling_generated}/probe_background.cpp"
        "${sampling_generated}/audited_math.cpp")
    target_include_directories(em_secondary_probe PRIVATE ../src/lib
        "${CMAKE_CURRENT_SOURCE_DIR}/sampling" "${MCRAY_NR_REFERENCE_DIR}")
    target_link_libraries(em_secondary_probe PRIVATE mcray)
    add_executable(secondary_oracle_tests sampling/secondary_oracle_tests.cpp)
    target_include_directories(secondary_oracle_tests PRIVATE ../src/lib)
    target_link_libraries(secondary_oracle_tests PRIVATE GTest::gtest_main mcray)
    gtest_discover_tests(secondary_oracle_tests TEST_PREFIX "unit.sampling."
        PROPERTIES TIMEOUT 30 LABELS "unit;sampling")
endif()
