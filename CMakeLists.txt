cmake_minimum_required(VERSION 3.0.0)
project(ellipsoidInclusion VERSION 0.1.0)

link_libraries(armadillo)
include(CTest)
enable_testing()
add_subdirectory(src)
add_subdirectory(test)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

add_test(NAME lb_contained
         COMMAND ellipsoidInclusionTest 1)

add_test(NAME lb_not_contained
         COMMAND ellipsoidInclusionTest 2)

add_test(NAME ells_contained
         COMMAND ellipsoidInclusionTest 3)

add_test(NAME ells_not_contained
         COMMAND ellipsoidInclusionTest 4)





# target_link_libraries(${PROJECT_TEST_NAME} gcov)
# setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
     