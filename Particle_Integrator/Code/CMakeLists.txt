cmake_minimum_required (VERSION 2.6)
project (ParticleIntegrator)

# The version number.
set (ParticleIntegrator_VERSION_MAJOR 0)
set (ParticleIntegrator_VERSION_MINOR 13)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/PIConfig.h.in"
  "${PROJECT_BINARY_DIR}/PIConfig.h"
  )

include_directories(${PROJECT_BINARY_DIR}) # PIConfig.h

include_directories(${PROJECT_SOURCE_DIR}/Inih)
include_directories(${PROJECT_SOURCE_DIR}/cspice/include)

IF(CMAKE_COMPILER_IS_GNUCC)
include_directories(/usr/include)
ENDIF()

add_executable(
	PIntegrator SourceCode.c
	resources/pintegrator.rc
)

option (WTIMING
    "Measure computing duration" ON)

if (WTIMING)
	add_definitions(-D__WTIMING)
	message("TIMING ACTIVE")
endif (WTIMING)

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
	SET(CUSTOM_C_FLAGS "/W4 /MD /O2 /Oi /Ot /Oy- /GL /Gy ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "/Wall /MDd /O2 /Oi /Ot /Oy- /GL /Gy ${CUSTOM_C_FLAG_OMP_DEBUG}")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CUSTOM_C_FLAGS "-O2 -funswitch-loops -Wall -Wpointer-arith -std=gnu99 ${CUSTOM_C_FLAG_OMP}")
	SET(CUSTOM_C_FLAGS_DEBUG "-Og -Wall -Wextra -Wpointer-arith -std=gnu99 -g ${CUSTOM_C_FLAG_OMP_DEBUG}")
ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CUSTOM_C_FLAGS}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CUSTOM_C_FLAGS_DEBUG}")
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${CUSTOM_C_FLAGS}")

message("> CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
message("> CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE}")
message("> CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG}")
message("> CMAKE_C_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO}")

set(CMAKE_BUILD_TYPE "Release")

#Full variable output for debugging
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
