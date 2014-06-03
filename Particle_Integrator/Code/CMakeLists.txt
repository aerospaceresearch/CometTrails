cmake_minimum_required (VERSION 2.6)
project (ParticleIntegrator)


# SOURCE DIR

set (PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR})


# VERSION

set (ParticleIntegrator_VERSION_MAJOR 0)
set (ParticleIntegrator_VERSION_MINOR 14)


# Passing on CMake Settings

configure_file (
	"${PROJECT_SOURCE_DIR}/PIConfig.h.in"
	"${PROJECT_BINARY_DIR}/PIConfig.h"
)


# EXECUTABLE

add_executable(
	PIntegrator SourceCode.c
	resources/pintegrator.rc
)


# OPTIONS

option (WTIMING
	"Measure computing duration" ON)
option (WSTEPINFO
	"Print time step extremes and number of steps" OFF)
option (ENDONTIME
	"Last step ends precisely on given end time" OFF)
option (PRD
	"Simulate Poynting-Robertson drag" ON)
option (SWD
	"Simulate solar wind drag" ON)
option (SaveRateOpt
	"Increase number of saved states when close to planets" ON)
if (MSVC)
	option (LinkStatic
		"Activate if cspice was compiled with /MT (pre-compiled versions are not)." OFF)
endif ()


# 32 or 64?

if ("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	set(CSPICE_FOLDER_NAME "cspice64")
elseif ("${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
	set(CSPICE_FOLDER_NAME "cspice32")
else ()
	message("Error: Unknown architecture.")
endif ()


# INCLUDE

include_directories(${PROJECT_BINARY_DIR}) # PIConfig.h

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/Inih)
include_directories(${PROJECT_SOURCE_DIR}/Integrators)
include_directories(${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/include)

if (CMAKE_COMPILER_IS_GNUCC)
include_directories(/usr/include)
endif ()


# DEFINE

if (WTIMING)
	add_definitions(-D__WTIMING)
	message("TIMING ACTIVE")
endif (WTIMING)
if (WSTEPINFO)
	add_definitions(-D__WTIMESTEP)
	message("TIME STEP EXTREMES OUTPUT ACTIVE")
endif (WSTEPINFO)
if (ENDONTIME)
	add_definitions(-D__ENDONTIME)
	message("ENDONTIME ACTIVE")
endif (ENDONTIME)
if (PRD)
	add_definitions(-D__PRD)
	message("PRD ACTIVE")
endif (PRD)
if (SWD)
	add_definitions(-D__SWD)
	message("SWD ACTIVE")
endif (SWD)
if (SaveRateOpt)
	add_definitions(-D__SaveRateOpt)
	message("SaveRateOpt ACTIVE")
endif (SaveRateOpt)


# LIBRARIES

if (MSVC)
	SET(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:LIBCMT;LIBCMTD)
	message("> CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
endif ()

add_subdirectory (${PROJECT_SOURCE_DIR}/Inih)
target_link_libraries (PIntegrator Inih)

if (MSVC)
find_library( CSPICE NAMES cspice
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib
)
find_library( CSUPPORT NAMES csupport
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib
)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT})
endif ()
if (CMAKE_COMPILER_IS_GNUCC)
find_library( CSPICE NAMES cspice.a
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib NO_DEFAULT_PATH 
	)
find_library( CSUPPORT NAMES csupport.a
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib NO_DEFAULT_PATH 
	)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT} m)
endif ()

target_link_libraries(PIntegrator ${LIBRARY_PATHS})


# FLAGS

if (MSVC)
	FIND_PACKAGE( OpenMP REQUIRED)
	if (OPENMP_FOUND)
		message("OPENMP FOUND")
		SET(CUSTOM_C_FLAG_OMP "/openmp")
		SET(CUSTOM_C_FLAG_OMP_DEBUG "/openmp")
	endif ()
endif ()
if (CMAKE_COMPILER_IS_GNUCC)
	message("GNUCC setting openmp flags")
	SET(CUSTOM_C_FLAG_OMP "-fopenmp")
	SET(CUSTOM_C_FLAG_OMP_DEBUG "-fopenmp")
endif ()

if (MSVC)
	SET(CUSTOM_C_FLAGS "/W4 /O2 /Oi /Ot /GL /Gy ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "/Wall /Oi /Ot /Oy- /GL /Gy ${CUSTOM_C_FLAG_OMP_DEBUG}")
	if (LinkStatic)
		message("Static linking (/MT) active")
		SET(CUSTOM_C_FLAGS "${CUSTOM_C_FLAGS} /MT")
		SET(CUSTOM_C_FLAGS_DEBUG "${CUSTOM_C_FLAGS_DEBUG} /MTd")
	elseif (LinkStatic)
		SET(CUSTOM_C_FLAGS "${CUSTOM_C_FLAGS} /MD")
		SET(CUSTOM_C_FLAGS_DEBUG "${CUSTOM_C_FLAGS_DEBUG} /MDd")
	endif (LinkStatic)
	SET(CUSTOM_C_FLAGS_RELWITHDEBINFO "${CUSTOM_C_FLAGS} /Oy-")
endif ()
if (CMAKE_COMPILER_IS_GNUCC)
	SET(CUSTOM_C_FLAGS "-O2 -funswitch-loops -Wall -Wpointer-arith -std=gnu99 -static-libgcc ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "-O0 -Wall -Wextra -Wpointer-arith -std=gnu99 -static-libgcc -ggdb3 ${CUSTOM_C_FLAG_OMP_DEBUG}")
endif ()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CUSTOM_C_FLAGS_DEBUG}")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CUSTOM_C_FLAGS_RELWITHDEBINFO}")

if (MSVC)
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /DRELTYPEREL")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /DRELTYPEDEB")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} /DRELTYPERWDI")
endif ()
if (CMAKE_COMPILER_IS_GNUCC)
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -D RELTYPEREL")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D RELTYPEDEB")
endif ()

message("> CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message("> CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
message("> CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")
if (MSVC)
message("> CMAKE_C_FLAGS_RELWITHDEBINFO: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
endif ()


# XP Platform toolset (edit manually)

IF(MSVC)
#	To set the platform toolset, uncomment the line with the desired toolkit, then delete the CMake build directory while CMake is closed.
#	Reset to default:
#		set(CMAKE_GENERATOR_TOOLSET "v120" CACHE STRING "Platform Toolset" FORCE) 
#	Win XP:
		set(CMAKE_GENERATOR_TOOLSET "v120_xp" CACHE STRING "Platform Toolset" FORCE) 
ENDIF(MSVC)
	message("> CMAKE_GENERATOR_TOOLSET: ${CMAKE_GENERATOR_TOOLSET}")


# DEBUG

#Full variable output for debugging
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
