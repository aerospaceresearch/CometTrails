cmake_minimum_required (VERSION 2.6)
project (ParticleIntegrator)


# SOURCE DIR

set (PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR})


# VERSION

set (ParticleIntegrator_VERSION_MAJOR 0)
set (ParticleIntegrator_VERSION_MINOR 13)

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


# INCLUDE

include_directories(${PROJECT_BINARY_DIR}) # PIConfig.h

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/Inih)
include_directories(${PROJECT_SOURCE_DIR}/Integrators)
include_directories(${PROJECT_SOURCE_DIR}/cspice/include)

IF(CMAKE_COMPILER_IS_GNUCC)
include_directories(/usr/include)
ENDIF()


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

IF(MSVC)
		SET(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:LIBCMT)
ENDIF()

add_subdirectory (${PROJECT_SOURCE_DIR}/Inih)
target_link_libraries (PIntegrator Inih)

IF(MSVC)
find_library( CSPICE NAMES cspice
			PATHS   ${PROJECT_SOURCE_DIR}/cspice/lib
)
find_library( CSUPPORT NAMES csupport
			PATHS   ${PROJECT_SOURCE_DIR}/cspice/lib
)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT})
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
find_library( CSPICE NAMES cspice.a
			PATHS   ${PROJECT_SOURCE_DIR}/cspice/lib NO_DEFAULT_PATH 
	)
find_library( CSUPPORT NAMES csupport.a
			PATHS   ${PROJECT_SOURCE_DIR}/cspice/lib NO_DEFAULT_PATH 
	)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT} m)
ENDIF()

target_link_libraries(PIntegrator ${LIBRARY_PATHS})


# FLAGS

IF(MSVC)
	FIND_PACKAGE( OpenMP REQUIRED)
	IF(OPENMP_FOUND)
		message("OPENMP FOUND")
		SET(CUSTOM_C_FLAG_OMP "/openmp")
		SET(CUSTOM_C_FLAG_OMP_DEBUG "/openmp")
	ENDIF()
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
	message("GNUCC setting openmp flags")
	SET(CUSTOM_C_FLAG_OMP "-fopenmp")
	SET(CUSTOM_C_FLAG_OMP_DEBUG "-fopenmp")
ENDIF()


IF(MSVC)
	SET(CUSTOM_C_FLAGS "/W4 /MD /O2 /Oi /Ot /GL /Gy ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "/Wall /MDd /Oi /Ot /Oy- /GL /Gy ${CUSTOM_C_FLAG_OMP_DEBUG}")
	SET(CUSTOM_C_FLAGS_RELWITHDEBINFO "${CUSTOM_C_FLAGS} /Oy-")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CUSTOM_C_FLAGS "-O2 -funswitch-loops -Wall -Wpointer-arith -std=gnu99 ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "-Og -Wall -Wextra -Wpointer-arith -std=gnu99 -g ${CUSTOM_C_FLAG_OMP_DEBUG}")
ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CUSTOM_C_FLAGS_DEBUG}")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CUSTOM_C_FLAGS_RELWITHDEBINFO}")

message("> CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message("> CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
message("> CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")
message("> CMAKE_C_FLAGS_RELWITHDEBINFO: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")


# DEBUG

#Full variable output for debugging
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
