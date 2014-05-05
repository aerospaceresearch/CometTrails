cmake_minimum_required (VERSION 2.6)
project (Inih)

add_library(Inih ini.c)

IF(MSVC)
	SET(PLATFORM_C_FLAGS "/W4 /MD /O2 /Oi /Ot /Oy- /GL /Gy")
	SET(PLATFORM_C_FLAGS_DEBUG "/Wall /MDd /O2 /Oi /Ot /Oy- /GL /Gy /Zi")
ENDIF()
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(PLATFORM_C_FLAGS "-O2 -funswitch-loops -Wall -Wpointer-arith -std=c99")
	SET(PLATFORM_C_FLAGS_DEBUG "-Og -Wall -Wextra -Wpointer-arith -std=c99 -g")
ENDIF()