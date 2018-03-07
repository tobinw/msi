#This is the top MSI CMakeList File for the Build

#Setting Version Number, Project Name
cmake_minimum_required (VERSION 3.4)
project (msi)
enable_language (Fortran)
set(MSI_VERSION 1.0.0)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#unless building shared libs, then select static libs
# if both static and shared libs are available
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wextra")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -Wall -Wextra")
endif()

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(ENABLE_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE string
    "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
    CACHE string
    "the command line flag to give process count to MPIRUN")

#Doxygen generation system
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
               ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

#Source and header files
set(SOURCES
  api/msi.cc
  src/msi_solver.cc
  src/msi_petsc.cc
  src/msi_trilinos.cc
)

set(HEADERS
  api/msi.h
  api/msi_types.h
  api/msi.f90
)

#  include/msi_las.h
#  include/msi_solver.h
#  include/msi_petsc.h
#  include/msi_trilinos.h
#  include/msi_field_storage.h
#  include/msi_field_op.h


set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

list(APPEND INCLUDES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                     $<INSTALL_INTERFACE:include>
                     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/api>)
set(TARGET_LIB_NAME msi)

find_package(Scorec REQUIRED)
list(APPEND INCLUDES ${SCOREC_INCLUDE_DIRS})
list(APPEND DEP_LIBS ${SCOREC_LIBRARIES})

if(ENABLE_PETSC)
  find_package(Hdf5 REQUIRED)
  find_package(PETSC REQUIRED)
  list(APPEND INCLUDES ${HDF5_INCLUDE_DIRS} ${PETSC_INCLUDE_DIRS})
  list(APPEND DEP_LIBS ${HDF5_LIBRARIES} ${PETSC_LIBRARIES})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMSI_PETSC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMSI_PETSC")
endif()

if(ENABLE_TRILINOS)
  find_package(Trilinos QUIET REQUIRED)
  find_package(Boost QUIET REQUIRED)
  list(APPEND INCLUDES ${TRILINOS_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS})
  list(APPEND DEP_LIBS ${TILINOS_LIBRARIES} ${BOOST_LIBRARIES})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMSI_TRILINOS")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMSI_TRILINOS")
  set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_trilinos)
else()
  set(ENABLE_PETSC "ON")
endif()

if(ENABLE_COMPLEX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPETSC_USE_COMPLEX")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPETSC_USE_COMPLEX")
  set(TARGET_LIB_NAME ${TARGET_LIB_NAME}_complex)
endif()

add_library(${TARGET_LIB_NAME} ${SOURCES})
target_include_directories(${TARGET_LIB_NAME} PUBLIC ${INCLUDES})
target_link_libraries(${TARGET_LIB_NAME} ${DEP_LIBS})

install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${TARGET_LIB_NAME}
        EXPORT ${TARGET_LIB_NAME}Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
export(EXPORT ${TARGET_LIB_NAME}Targets
       FILE "${CMAKE_BINARY_DIR}/msi/${TARGET_LIB_NAME}Targets.cmake")
set(ConfigPackageLocation lib/cmake/msi)
install(EXPORT ${TARGET_LIB_NAME}Targets
        DESTINATION ${ConfigPackageLocation})
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/msi/${TARGET_LIB_NAME}ConfigVersion.cmake"
                                 VERSION ${MSI_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/msiConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/msi/${TARGET_LIB_NAME}Config.cmake" @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/msi/${TARGET_LIB_NAME}Config.cmake ${CMAKE_CURRENT_BINARY_DIR}/msi/${TARGET_LIB_NAME}ConfigVersion.cmake
        DESTINATION ${ConfigPackageLocation}
        COMPONENT Devel)

if (ENABLE_TESTING)
  if (ENABLE_PETSC)
  add_executable(test_petsc test/petsc/main.cc)
  target_link_libraries(test_petsc ${TARGET_LIB_NAME} ${DEP_LIBS})
  message(STATUS "test_petsc link libraries: ${TARGET_LIB_NAME} ${DEP_LIBS}")
  INSTALL(TARGETS test_petsc RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()
if (ENABLE_TRILINOS)
    add_executable(epetra test/epetra/main.cc)
    target_link_libraries(epetra ${DEP_LIBS})
    target_link_libraries(epetra ${TARGET_LIB_NAME})
    INSTALL(TARGETS epetra RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endif()
endif() #enable_testing

#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)

message("ENV_FLAGS = $ENV{FFLAGS}")
message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_Fortran_COMPILER_INIT = ${CMAKE_Fortran_COMPILER_INIT}")
message("CMAKE_Fortran_COMPILER_FULLPATH = ${CMAKE_Fortran_COMPILER_FULLPATH}")
message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
message("CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")


