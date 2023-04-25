cmake_minimum_required( VERSION 3.16 )

enable_language( Fortran )

option( USE_MPI    OFF )
option( USE_OPENMP OFF )

set( PROJECT_TARGET odeExercise )
project( ${PROJECT_TARGET} )
set( CMAKE_BUILD_TYPE Debug )

# Custom netCDF finder
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake )
find_package( netCDF REQUIRED ) # NO_DEFAULT_PATH PATHS ${PROJECT_SOURCE_DIR} )

# Create exec
add_executable( ${PROJECT_TARGET}
                solver.f90
                equations.f90
                observer.f90
                )


if ( USE_MPI )
  find_package( MPI REQUIRED COMPONENTS Fortran )
  target_compile_definitions( ${PROJECT_TARGET} PRIVATE -DUSE_MPI=1 )
endif()

if ( USE_OPENMP )
  find_package( OpenMP REQUIRED COMPONENTS Fortran )
  if ( NOT OpenMP_Fortran_FOUND )
    message( "OpenMP_Fortran_FOUND => ${OpenMP_Fortran_FOUND}" )
  endif()
endif()

# https://stackoverflow.com/a/53155812
set( Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID} )
message( STATUS "Set Fortran_COMPILER_ID to : ${Fortran_COMPILER_ID}" )
target_compile_options( ${PROJECT_TARGET} 
                        PRIVATE
                          -cpp
                          # $<$<STREQUAL:"${Fortran_COMPILER_ID}","GNU">:-fcheck=all>
                          # $<$<STREQUAL:"${Fortran_COMPILER_ID}","Intel">:-check all>
                        )

target_link_libraries(  ${PROJECT_TARGET}
                          PRIVATE
                            ${netCDF_LIBRARIES}
                            $<$<BOOL:${USE_MPI}>:$<TARGET_NAME_IF_EXISTS:MPI::MPI_Fortran>>
                            $<$<BOOL:${USE_OPENMP}>:$<TARGET_NAME_IF_EXISTS:OpenMP::OpenMP_Fortran>>
                        )

target_include_directories( ${PROJECT_TARGET}
                            PRIVATE
                              ${netCDF_INCLUDE_DIRS}
                            )

# message( "netCDF_INCLUDE_DIRS => ${netCDF_INCLUDE_DIRS}" )
# message( "netCDF_LIBRARIES    => ${netCDF_LIBRARIES}" )

# message( "netCDF_INCLUDE_DIR  => ${netCDF_INCLUDE_DIR}" )
# message( "netCDF_CLIBS        => ${netCDF_CLIBS}" )
# message( "netCDF_CXXLIBS      => ${netCDF_CXXLIBS}" )
# message( "netCDF_FLIBS        => ${netCDF_FLIBS}" )

install(
        TARGETS odeExercise
        RUNTIME DESTINATION bin/
        ARCHIVE DESTINATION lib/
        LIBRARY DESTINATION lib/
        )