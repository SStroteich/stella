# FindSLEPc.cmake
# Find the SLEPc library

# Ensure SLEPC_DIR and PETSC_ARCH are set
if(NOT DEFINED ENV{SLEPC_DIR})
  message(FATAL_ERROR "SLEPC_DIR environment variable is not set")
endif()
if(NOT DEFINED ENV{PETSC_ARCH})
  message(FATAL_ERROR "PETSC_ARCH environment variable is not set")
endif()

# Set the paths
set(SLEPC_DIR $ENV{SLEPC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})
set(SLEPC_INCLUDE_DIR "${SLEPC_DIR}/include")
set(SLEPC_ARCH_INCLUDE_DIR "${SLEPC_DIR}/${PETSC_ARCH}/include")
set(SLEPC_LIB_DIR "${SLEPC_DIR}/${PETSC_ARCH}/lib")

# Locate SLEPc headers
find_path(SLEPC_INCLUDES slepc.h
  HINTS ${SLEPC_INCLUDE_DIR} ${SLEPC_ARCH_INCLUDE_DIR}
)

# Locate SLEPc libraries
find_library(SLEPC_LIBRARIES NAMES slepc
  HINTS ${SLEPC_LIB_DIR}
)

# Handle the standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc DEFAULT_MSG SLEPC_LIBRARIES SLEPC_INCLUDES)

# Set the include directories and libraries
if(SLEPC_FOUND)
  set(SLEPC_INCLUDE_DIRS ${SLEPC_INCLUDE_DIR} ${SLEPC_ARCH_INCLUDE_DIR})
  set(SLEPC_LIBRARIES ${SLEPC_LIBRARIES})
endif()

