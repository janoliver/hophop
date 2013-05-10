# This module defines
# GSL_INCLUDE_DIR, where to find include files, etc.
# GSL_LIBRARIES, the libraries to link against.
# GSL_STATIC_LIBRARY_PATH
# GSL_FOUND.

# also defined, but not for general use are
# GSL_LIBRARY, where to find the CUnit library.

message(STATUS "Searching for gsl library")

find_path(GSL_INCLUDE_DIR gsl)

find_library(GSL_LIBRARY gsl)
find_library(GSLCBLAS_LIBRARY gslcblas)

if(GSL_INCLUDE_DIR)
  if(GSL_LIBRARY)
    set(GSL_FOUND TRUE)
    set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
  endif(GSL_LIBRARY)
endif(GSL_INCLUDE_DIR)

if (GSL_FOUND)
   if (NOT GSL_FIND_QUIETLY)
      message(STATUS "Found GSL: ${GSL_LIBRARIES}")
   endif (NOT GSL_FIND_QUIETLY)
else (GSL_FOUND)
   if (GSL_FIND_REQUIRED)
      message(SEND_ERROR "Could NOT find GSL")
   endif (GSL_FIND_REQUIRED)
endif (GSL_FOUND)

mark_as_advanced (
  GSL_INCLUDE_DIR 
  GSL_LIBRARIES
)
