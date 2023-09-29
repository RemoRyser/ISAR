cmake_minimum_required (VERSION 2.6)
project (gridSimulationAFKtest)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)


# enable fortran, c, and c++ language
enable_language(C CXX)

##define compiler version
#check cxx standard
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()


###################### gsl include #####################
include_directories( /usr/local/gsl/1.16-3/include )
find_package( PkgConfig REQUIRED)
pkg_check_modules( gsl REQUIRED gsl )
set(ENV{PATH} "/usr/local/gsl/1.16-3/lib/")

#########################################################

###############sundials###################
#For the shared library:
set ( PROJECT_LINK_LIBS sundials_cvode sundials_nvecserial)
link_directories( /data/idiv_brose/deps/lib/ )
include_directories( /data/idiv_brose/deps/include/)

##########################################
#include_directories(*.h)

add_executable(simulation main.cpp)
target_link_libraries(simulation ${PROJECT_LINK_LIBS} )
target_link_libraries( simulation ${gsl_LIBRARIES} )
