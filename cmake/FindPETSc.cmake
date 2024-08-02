# FindPETSc.cmake
# Find the PETSc library

# Ensure PETSC_DIR and PETSC_ARCH are set
if(NOT DEFINED ENV{PETSC_DIR})
  message(FATAL_ERROR "PETSC_DIR environment variable is not set")
endif()
if(NOT DEFINED ENV{PETSC_ARCH})
  message(FATAL_ERROR "PETSC_ARCH environment variable is not set")
endif()

# Set the paths
set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})
set(PETSC_INCLUDE_DIR "${PETSC_DIR}/include")
set(PETSC_ARCH_INCLUDE_DIR "${PETSC_DIR}/${PETSC_ARCH}/include")
set(PETSC_LIB_DIR "${PETSC_DIR}/${PETSC_ARCH}/lib")

# Locate PETSc headers
find_path(PETSC_INCLUDES petsc.h
  HINTS ${PETSC_INCLUDE_DIR} ${PETSC_ARCH_INCLUDE_DIR}
)

# Locate PETSc libraries
find_library(PETSC_LIBRARIES NAMES petsc
  HINTS ${PETSC_LIB_DIR}
)

# Handle the standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDES)

# Set the include directories and libraries
if(PETSC_FOUND)
  set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} ${PETSC_ARCH_INCLUDE_DIR})
  set(PETSC_LIBRARIES ${PETSC_LIBRARIES})
endif()

