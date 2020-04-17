cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

#set(CMAKE_VERBOSE_MAKEFILE ON) 

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH}) 
project(SkyNetDriver CXX C)

# Add the executable
add_executable(postprocess src/postprocess.cpp)
add_executable(fulltraj src/fulltraj.cpp)
add_executable(reduce_reactions src/reduce_reactions.cpp)


# Required packages
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

 
find_package(SkyNet REQUIRED)



include_directories(${SKYNET_INCLUDE_DIRS})
target_link_libraries(postprocess ${SKYNET_LIBRARIES}) 
target_link_libraries(reduce_reactions ${SKYNET_LIBRARIES}) 
target_link_libraries(fulltraj ${MPI_LIBRARIES})
target_link_libraries(fulltraj ${SKYNET_LIBRARIES}) 

#find_package(PythonLibs 2.6 REQUIRED)
#include_directories(${PYTHON_INCLUDE_DIRS})
#target_link_libraries(rprocess ${PythonLibs_LIBRARIES}) 
 
#find_package(BLAS REQUIRED) 
#target_link_libraries(rprocess ${BLAS_LIBRARIES}) 
#find_package(LAPACK REQUIRED) 
#target_link_libraries(rprocess ${LAPACK_LIBRARIES}) 

#find_package(Boost REQUIRED system filesystem serialization)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(rprocess ${Boost_LIBRARIES}) 


# enable warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

# enable C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Enable floating point exception handling 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fnon-call-exceptions")

# generate debugging symbols for release and debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# the compiler flags
get_directory_property(OUT_VAR LINK_DIRECTORIES) 
message(STATUS "")
message(STATUS "")
message(STATUS "Configuration options:")
message(STATUS "")
message(STATUS "              CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "         CMAKE_CXX_LINK_FLAGS: " ${CMAKE_CXX_LINK_FLAGS})
message(STATUS "        CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "      CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "             LINK_DIRECTORIES: " ${OUT_VAR})
message(STATUS "             BOOST_DIRECTORIES: " ${Boost_LIBRARIES})

message(STATUS "")
message(STATUS "")

# Add subdirectories
# add_subdirectory(tests)


