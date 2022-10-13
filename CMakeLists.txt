cmake_minimum_required(VERSION 3.0.0)
project(ellipsoidInclusion VERSION 0.1.0)

link_libraries(armadillo)
include(CTest)
enable_testing()

# file (GLOB_RECURSE Foo_SOURCES CONFIGURE_DEPENDS "ellincheck.cpp")
# file (GLOB_RECURSE Foo_HEADERS CONFIGURE_DEPENDS "ellincheck.hpp")
add_executable(ellipsoidInclusion ellincheck.cpp ellincheck.hpp main.cpp)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
