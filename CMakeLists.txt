cmake_minimum_required (VERSION 2.6)
project(LongestPath CXX)
set(LongestPath_VERSION_MAJOR 1)
set(LongestPath_VERSION_MINOR 0)

find_library(LIBGMP gmp REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(${Boost_INCLUDE_DIRS})
add_executable(longest_path src/main.cpp src/parse_graph.cpp)
set_target_properties(longest_path PROPERTIES COMPILE_FLAGS "-std=c++11 -Wall -pedantic -O3")
target_link_libraries(longest_path ${Boost_LIBRARIES} ${LIBGMP})
install(TARGETS longest_path DESTINATION bin)

# Testing

enable_testing()

macro(do_test arg)
  add_test(test_${arg} longest_path --input-file ${PROJECT_SOURCE_DIR}/tests/${arg}.input 2>/dev/null | diff - ${PROJECT_SOURCE_DIR}/tests/${arg}.output)
endmacro(do_test)

do_test(3x3_sq)
do_test(3x4_sq)
do_test(3x5_sq)
do_test(3x6_sq)
do_test(3x7_sq)
do_test(3x8_sq)
do_test(4x4_sq)
do_test(4x5_sq)
do_test(4x6_sq)
do_test(4x7_sq)
do_test(4x8_sq)
do_test(5x5_sq)
do_test(5x6_sq)
do_test(5x7_sq)
do_test(5x8_sq)
do_test(6x6_sq)
do_test(6x7_sq)
do_test(6x8_sq)

# http://stackoverflow.com/questions/733475/cmake-ctest-make-test-doesnt-build-tests
# http://public.kitware.com/Bug/view.php?id=8774
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS longest_path)

# Packaging

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Computes the Longest Path Polynomial")
set(CPACK_PACKAGE_CONTACT "Andrea Bedini <andrea@andreabedini.com>")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${LongestPath_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${LongestPath_VERSION_MINOR}")
include(CPack)
