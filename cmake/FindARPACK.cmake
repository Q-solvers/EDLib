SET(ARPACK_SEARCH_PATHS ${ARPACK_DIR})

IF (DEFINED ARPACK_SEARCH_PATHS)
    FIND_LIBRARY(ARPACK_LIB NAMES arpack HINTS ${ARPACK_SEARCH_PATHS} PATH_SUFFIXES lib lib64)
ELSE (DEFINED ARPACK_SEARCH_PATHS)
    FIND_LIBRARY(ARPACK_LIB NAMES arpack HINTS $ENV{ARPACK_DIR} PATH_SUFFIXES lib lib64)
ENDIF (DEFINED ARPACK_SEARCH_PATHS)

if(USE_MPI)
  IF (DEFINED ARPACK_SEARCH_PATHS)
    FIND_LIBRARY(PARPACK_LIB NAMES parpack HINTS ${ARPACK_SEARCH_PATHS} PATH_SUFFIXES lib lib64)
  ELSE (DEFINED ARPACK_SEARCH_PATHS)
    FIND_LIBRARY(PARPACK_LIB NAMES parpack HINTS $ENV{ARPACK_DIR} PATH_SUFFIXES lib lib64)
  ENDIF (DEFINED ARPACK_SEARCH_PATHS)
endif (USE_MPI)

SET(ARPACK_FOUND FALSE)
IF (ARPACK_LIB)
    SET(ARPACK_FOUND TRUE)
    IF (PARPACK_LIB)
        SET(PARPACK_FOUND TRUE)
        MARK_AS_ADVANCED(PARPACK_LIB)
    ENDIF (PARPACK_LIB)
    MARK_AS_ADVANCED(ARPACK_LIB)
ENDIF (ARPACK_LIB)

IF (ARPACK_FOUND)
    IF (NOT ARPACK_LIB_FIND_QUIETLY)
        MESSAGE(STATUS "Found Arpack : ${ARPACK_LIB}")
        GET_FILENAME_COMPONENT(ARPACK_PATH ${ARPACK_LIB} PATH CACHE)
        SET(ARPACK_INCLUDE_DIR ${ARPACK_PATH}/../include CACHE FILEPATH "ARPACK include directory.")
    ENDIF (NOT ARPACK_LIB_FIND_QUIETLY)
ELSE(ARPACK_FOUND)
    IF (ARPACK_LIB_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find Arpack")
    ENDIF (ARPACK_LIB_FIND_REQUIRED)
ENDIF (ARPACK_FOUND)
