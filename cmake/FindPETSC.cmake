find_package(PkgConfig REQUIRED QUIET)

list(APPEND CMAKE_PREFIX_PATH "${PETSC_DIR}/${PETSC_ARCH}")

if(PETSC_FIND_REQUIRED)
  set(_PETSC_OPTS "REQUIRED")
endif()
if(PETSC_FIND_QUIETLY)
  set(_PETSC_OPTS "QUIET")
endif()
if(PETSC_FIND_REQUIRED AND PETSC_FIND_QUIETLY)
  set(_PETSC_OPTS "REQUIRED QUIET")
endif()

pkg_check_modules(PETSC ${_PETSC_OPTS} IMPORTED_TARGET PETSc)

if(PETSC_FOUND)
  set(PETSC_INCLUDE_DIRS "${PETSC_STATIC_INCLUDE_DIRS}")
  set(PETSC_LIBRARIES "${PETSC_STATIC_LDFLAGS}")
  if("-lzoltan" IN_LIST PETSC_LIBRARIES)
    message(STATUS "Found Zoltan library as PETSc dependency.")
    set(ZOLTAN_FOUND true)
    set(ZOLTAN_LIBRARIES )
    set(ZOLTAN_INCLUDE_DIRS )
  endif()
  if("-lparmetis" IN_LIST PETSC_LIBRARIES)
    message(STATUS "Found Parmetis library as PETSc dependency.")
    set(PARMETIS_FOUND true)
    set(PARMETIS_LIBRARIES )
    set(PARMETIS_INCLUDE_DIRS )
  endif()
  if("-lmetis" IN_LIST PETSC_LIBRARIES)
    message(STATUS "Found Metis library as PETSc dependency.")
    set(METIS_FOUND true)
    set(METIS_LIBRARIES )
    set(METIS_INCLUDE_DIRS )
  endif()
endif(PETSC_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC DEFAULT_MSG PETSC_LIBRARIES PETSC_INCLUDE_DIRS)
mark_as_advanced(PETSC_LIBRARIES PETSC_INCLUDE_DIRS)