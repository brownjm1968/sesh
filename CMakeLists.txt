cmake_minimum_required(VERSION 3.5)

project(SESH)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form  -fimplicit-none") #-std=f2008
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
file(GLOB_RECURSE sources  src/mc-routines.f90
                           src/physics.f90
                           src/sesh.f90 )
file(GLOB_RECURSE testsrc  test/test-sesh.f90 )
add_executable(sesh ${sources})
add_executable(test-sesh ${testsrc} src/physics.f90 )

# target_link_libraries(test-sesh sources)
