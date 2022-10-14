cmake_minimum_required(VERSION 3.0.0)
project(ellipsoidInclusion VERSION 0.1.0)

link_libraries(armadillo)
include(CTest)
enable_testing()
add_subdirectory(ellincheck)


add_executable(ellipsoidInclusion main.cpp)
target_link_libraries(ellipsoidInclusion PRIVATE ellincheck)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

add_test(NAME ellipsoidInclusion
         COMMAND ellipsoidInclusion)



set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
include(CodeCoverage)
append_coverage_compiler_flags()

# target_link_libraries(${PROJECT_TEST_NAME} gcov)
# setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
     