cmake_minimum_required(VERSION 2.6)
    
project(aocommon)

add_executable(uvector main.cpp)
set_target_properties(uvector PROPERTIES COMPILE_FLAGS "-std=c++11 -Wall -pedantic -O3 -march=native")

add_executable(tcompare tcompare.cpp stopwatch.cpp)
set_target_properties(tcompare PROPERTIES COMPILE_FLAGS "-std=c++11 -Wall -pedantic -O3 -march=native")

add_executable(tlane03 tlane03.cpp stopwatch.cpp)
set_target_properties(tlane03 PROPERTIES COMPILE_FLAGS "-std=c++11 -Wall -pedantic -O3 -march=native")
target_link_libraries(tlane03 boost_thread boost_system)

add_executable(tlane11 tlane11.cpp stopwatch.cpp)
set_target_properties(tlane11 PROPERTIES COMPILE_FLAGS "-std=c++11 -Wall -pedantic -O3 -march=native -pthread")
set_target_properties(tlane11 PROPERTIES LINK_FLAGS "-pthread")

install(TARGETS uvector RUNTIME DESTINATION bin)
install(FILES lane.h lane_03.h lane_11.h uvector.h DESTINATION include)

find_package(Doxygen)

if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
	${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM)

endif(DOXYGEN_FOUND)
