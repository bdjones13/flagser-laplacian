cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(flagser-laplacian)

enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

find_package(Threads REQUIRED)

# ##########################################################################
# ######    CHANGE THIS DEPENDING ON YOUR MACHINE   ######################## 
# set(Matlab_ROOT_DIR "/usr/local/MATLAB/R2022b")###########################
# ##########################################################################

# find_package(Matlab)
# find_library(LIB_MATLAB_DATA_ARRAY 
# 	NAMES MatlabDataArray
# 	PATHS ${Matlab_ROOT_DIR}/extern/bin/glnxa64
# 	NO_DEFAULT_PATH
#   REQUIRED)

# find_library(LIB_MATLAB_ENGINE
# 	NAMES MatlabEngine
# 	PATHS ${Matlab_ROOT_DIR}/extern/bin/glnxa64
#         NO_DEFAULT_PATH
#         REQUIRED)

# include_directories(${Matlab_INCLUDE_DIRS})


IF(MSVC)
  set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_FLAGS} -Ox -D NDEBUG")
ELSE()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pthread") # removed -pedantic-errors because Matlab throws some
  set(CMAKE_CXX_DEBUG_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fno-omit-frame-pointer -fsanitize=memory")
  set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -D NDEBUG")
ENDIF()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_POSITION_INDEPENDANT}")
set(CMAKE_VS_PLATFORM_NAME "x64")
set(CMAKE_VS_PLATFORM_TOOLSET_HOST_ARCHITECTURE "x64")


IF(NOT MSVC)
add_executable(flagser-laplacian src/flagser-laplacian.cpp)
target_compile_definitions(flagser-laplacian PRIVATE USE_COEFFICIENTS=1 REAL_COEFFICIENTS=1)
set_target_properties(flagser-laplacian PROPERTIES LINKER_LANGUAGE CXX)

ELSE()
  message("On Windows, compilation with coefficients is currently not supported.")
ENDIF()

message("")
message("** CMAKE_BUILD_TYPE: \t\t\t${CMAKE_BUILD_TYPE}")
message("** CMAKE_CXX_FLAGS: \t\t\t${CMAKE_CXX_FLAGS}")
message("** CMAKE_CXX_DEBUG_FLAGS: \t\t\t${CMAKE_CXX_DEBUG_FLAGS}")
message("** CMAKE_CXX_RELEASE_FLAGS: \t\t\t${CMAKE_CXX_RELEASE_FLAGS}")
message("** CMAKE_CXX_LINK_FLAGS: \t\t${CMAKE_CXX_LINK_FLAGS}")
message("")

find_library(PERSISTENTLAPLACIANS_LIB PersistentLaplacians HINTS /usr/local/lib /usr/local/include /mnt/home/jones657/customlib/include /mnt/home/jones657/customlib/lib)
message(PERSISTENTLAPLACIANS_LIB="${PERSISTENTLAPLACIANS_LIB}")

find_package(Eigen3 3.4 REQUIRED HINTS /usr/local/lib /usr/local/include /opt/software/Eigen/3.4.0-GCCcore-11.3.0/include)
link_directories(${Eigen_INCLUDE_DIRS})

include_directories(/usr/local/include /mnt/home/jones657/customlib/include/ ${Eigen_INCLUDE_DIRS} /opt/software/Eigen/3.4.0-GCCcore-11.3.0/include/ )

# add_executable(flagser-laplacian src/flagser-laplacian.cpp)
target_compile_definitions(flagser-laplacian PRIVATE USE_COEFFICIENTS=1 REAL_COEFFICIENTS=1)
set_target_properties(flagser-laplacian PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(flagser-laplacian ${HDF5_LIBRARIES} ${Eigen_LIBRARIES} ${PERSISTENTLAPLACIANS_LIB})


IF(NOT MSVC)
  install(TARGETS flagser-laplacian RUNTIME DESTINATION bin) 
# ELSE()
  # install(TARGETS flagser flagser-count RUNTIME DESTINATION bin)
ENDIF()

IF(NOT MSVC)
add_custom_command( # Ben Jones 2023-03-01
        TARGET flagser-laplacian POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:flagser-laplacian> ${CMAKE_SOURCE_DIR})

ENDIF()

