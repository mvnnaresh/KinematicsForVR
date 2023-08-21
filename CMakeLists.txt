cmake_minimum_required(VERSION 2.8.12)
set(QT_MINIMUM_VERSION 5.9.0)

project(VR_Kinematics)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

message("CMAKE_PROJECT_NAME = ${CMAKE_PROJECT_NAME}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

#############################################
#        		DEPENDENCIES			    #
#############################################

option(WITH_VISP "Build with VISP support" TRUE)

#Preprocessor flags
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

find_package(VISP REQUIRED)
find_package(orocos_kdl REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

file(GLOB SOURCES
        "${PROJECT_SOURCE_DIR}/include/*.h"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.c"
        )

file(GLOB HEADERS
	"${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
)

add_executable(${PROJECT_NAME} ${SOURCES})
#add_library(${CMAKE_PROJECT_NAME} SHARED ${SOURCES})

message(STATUS "Eigen INCLUDE_DIRECTORIES $ENV{EIGEN3_INCLUDE_DIR}")

target_compile_definitions(
		${CMAKE_PROJECT_NAME}
		PUBLIC
)

target_include_directories(
		${CMAKE_PROJECT_NAME}
		PUBLIC
		${VISP_INCLUDE_DIRS}
		${orocos_kdl_INCLUDE_DIRS}
	)
message("Orocos include: ${orocos_kdl_INCLUDE_DIRS}")

	target_link_libraries(
		${CMAKE_PROJECT_NAME}
		${VISP_LIBRARIES}	
		${orocos_kdl_LIBRARIES}	)

# note that it is not CMAKE_INSTALL_PREFIX we are checking here
set(CMAKE_INSTALL_PREFIX "C:/Users/Naresh-Admin/Desktop/VR_KIN")
if(DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX is not set\n"
        "Default value: ${CMAKE_INSTALL_PREFIX}\n"
        "Will set it to ${CMAKE_SOURCE_DIR}/install"
    )
    set(CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/install"
        CACHE PATH "Where the library will be installed to" FORCE
    )
else()
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX was already set\n"
        "Current value: ${CMAKE_INSTALL_PREFIX}"
    )
endif()

# without it public headers won't get installed
set(public_headers
    "${PROJECT_SOURCE_DIR}/include/*.h"
)
# note that ${public_headers} has to be in quotes
set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${HEADERS})

install(TARGETS ${CMAKE_PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
