cmake_minimum_required(VERSION 2.8)

#----------------------------------------------------------------
#-- may have here to precise the path to gcc and g++ if necessary
#----------------------------------------------------------------
set(CMAKE_C_COMPILER   gcc)
set(CMAKE_CXX_COMPILER g++)



project(hj_Distrib)
  
message ( "-- OS: ${CMAKE_SYSTEM}")

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(bitness 64)
  set(BITFLAG  "")
  message(".. bits are ${bitness} ")
else()
  set(bitness 32)
  set(BITFLAG  " -m32 ")
  message(".. bits are ${bitness} ")
endif()


## To enforce compilation with 32 bits: uncomment the following two ligns
#set(BITFLAG  " -m32 ")
#message(" enforce BITFLAG = ${BITFLAG}")

 
SET(SRC "${CMAKE_SOURCE_DIR}/src")

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

#if(${bitness} EQUAL 64)
#  set(CMAKE_C_FLAGS "   -O3 -ffast-math -march=opteron -funroll-loops")
#  set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})
#  link_directories(${CMAKE_SOURCE_DIR}/lib64)
#endif()

#if(${bitness} EQUAL 32)
set(CMAKE_C_FLAGS "${BITFLAG} -O3 -fopenmp")
message("option:  ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})
link_directories(${CMAKE_SOURCE_DIR}/lib)
#endif()
 
FILE(GLOB_RECURSE source_files ${CMAKE_SOURCE_DIR}/src/*.cpp) # TO COMPILE ALL .CPP FILES of src/
message("src du projet: ${source_files}")
ADD_EXECUTABLE(exe ${source_files})
#if(${bitness} EQUAL 64)
#  TARGET_LINK_LIBRARIES(exe hj64)
#endif()
#if(${bitness} EQUAL 32)
TARGET_LINK_LIBRARIES(exe hj)
#endif()

