cmake_minimum_required (VERSION 2.6)
project (ParticleAnalyser)


# SOURCE DIR

set (PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR})


# VERSION

set (ParticleAnalyser_VERSION_MAJOR 0)
set (ParticleAnalyser_VERSION_MINOR 01)

# Passing on CMake Settings

configure_file (
	"${PROJECT_SOURCE_DIR}/PAConfig.h.in"
	"${PROJECT_BINARY_DIR}/PAConfig.h"
)


# EXECUTABLE

add_executable(
	PAnalyser SourceCode.c
	resources/panalyser.rc
)


# OPTIONS

option (WTIMING
	"Measure computing duration" ON)
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

include_directories(${PROJECT_BINARY_DIR}) # PAConfig.h

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/include)

IF(CMAKE_COMPILER_IS_GNUCC)
include_directories(/usr/include)
ENDIF()


# DEFINE

if (WTIMING)
	add_definitions(-D__WTIMING)
	message("TIMING ACTIVE")
endif (WTIMING)


# LIBRARIES

IF(MSVC)
	SET(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:LIBCMT;LIBCMTD)
	message("> CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
ENDIF()

IF(MSVC)
find_library( CSPICE NAMES cspice
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib
)
find_library( CSUPPORT NAMES csupport
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib
)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT})
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
find_library( CSPICE NAMES cspice.a
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib NO_DEFAULT_PATH 
	)
find_library( CSUPPORT NAMES csupport.a
			PATHS   ${PROJECT_SOURCE_DIR}/${CSPICE_FOLDER_NAME}/lib NO_DEFAULT_PATH 
	)
SET(LIBRARY_PATHS ${CSPICE};${CSUPPORT} m)
ENDIF()

target_link_libraries(PAnalyser ${LIBRARY_PATHS})


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
	SET(CUSTOM_C_FLAGS "/W4 /O2 /Oi /Ot /GL /Gy ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "/Wall /Oi /Ot /Oy- /GL /Gy ${CUSTOM_C_FLAG_OMP_DEBUG}")
	if (LinkStatic)
		SET(CUSTOM_C_FLAGS "${CUSTOM_C_FLAGS} /MT")
		SET(CUSTOM_C_FLAGS_DEBUG "${CUSTOM_C_FLAGS_DEBUG} /MTd")
	elseif (LinkStatic)
		SET(CUSTOM_C_FLAGS "${CUSTOM_C_FLAGS} /MD")
		SET(CUSTOM_C_FLAGS_DEBUG "${CUSTOM_C_FLAGS_DEBUG} /MDd")
	endif (LinkStatic)
	SET(CUSTOM_C_FLAGS_RELWITHDEBINFO "${CUSTOM_C_FLAGS} /Oy-")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CUSTOM_C_FLAGS "-O2 -funswitch-loops -Wall -Wpointer-arith -std=gnu99 ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "-O0 -Wall -Wextra -Wpointer-arith -std=gnu99 -g ${CUSTOM_C_FLAG_OMP_DEBUG}")
ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CUSTOM_C_FLAGS_DEBUG}")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CUSTOM_C_FLAGS_RELWITHDEBINFO}")

message("> CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message("> CMAKE_C_FLAGS_RELEASE: ${CMAKE_C_FLAGS_RELEASE}")
message("> CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")
message("> CMAKE_C_FLAGS_RELWITHDEBINFO: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")


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
