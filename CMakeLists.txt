cmake_minimum_required(VERSION 2.8)
 
project(Demons)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

ADD_EXECUTABLE(DemonsRegistration DemonsRegistration-itk3.cxx)
TARGET_LINK_LIBRARIES( DemonsRegistration ${ITK_LIBRARIES})

ADD_EXECUTABLE(DemonsWarp DemonsWarp.cxx)
TARGET_LINK_LIBRARIES( DemonsWarp ${ITK_LIBRARIES})


