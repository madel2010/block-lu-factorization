PROJECT(blockklu)

cmake_minimum_required(VERSION 2.6)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

SET(CMAKE_CXX_FLAGS "-O0 -g -fopenmp -Wall")
SET(CMAKE_C_FLAGS "-O0 -g -fopenmp -Wall")
SET(LIBS gomp)

SET(INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/Include 
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/AMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/BTF/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CCOLAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/COLAMD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CSparse/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/CXSparse/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/LDL/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/UMFPACK/Include
		  ${PROJECT_SOURCE_DIR}/SuiteSparse/SuiteSparse_config)


set(LIBSPARSE_LIBS
	"${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/Lib/libklu.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/UMFPACK/Lib/libumfpack.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/AMD/Lib/libamd.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Lib/libcholmod.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/COLAMD/Lib/libcolamd.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/CCOLAMD/Lib/libccolamd.a"
	"${PROJECT_SOURCE_DIR}/SuiteSparse/CAMD/Lib/libcamd.a"
        "${PROJECT_SOURCE_DIR}/SuiteSparse/BTF/Lib/libbtf.a" 
)

#######START: Find Lapack #####################
find_package(PkgConfig)

find_library(HAVE_LAPACK lapack)
if(HAVE_LAPACK)
    SET(LIBS ${LIBS} lapack)
else()
    ##TRY USING PKGCONFIG
    pkg_check_modules(LAPACK lapack)
    if(LAPACK_FOUND)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${LAPACK_CFLAGS}")
        set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${LAPACK_CFLAGS}")
        set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
        include_directories(${LAPACK_INCLUDEDIR})
        link_directories(${LAPACK_LIBDIR})
    else()
	unset(LAPACK_FOUND CACHE)
	message(FATAL_ERROR "Can not find lapack library:")
    endif()
	
endif()

find_library(HAVE_BLAS blas)
if(HAVE_BLAS)
    SET(LIBS ${LIBS} blas)
else()
    ##try using pckg_config
    pkg_check_modules(BLAS blas)
    if(BLAS_FOUND)
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${BLAS_CFLAGS}")
        set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${BLAS_CFLAGS}")
        set(LIBS ${LIBS} ${BLAS_LIBRARIES})
        include_directories(${BLAS_INCLUDEDIR})
        link_directories(${BLAS_LIBDIR})
    else()
	unset(LAPACK_FOUND CACHE)
	message(FATAL_ERROR "Can not find blas library:")
    endif()
    
endif()
######END: Find Lapack #####################


#add definitions, compiler switches, etc.
ADD_DEFINITIONS(-Wall -O0 -g)

#ADD_SUBDIRECTORY(src)

include_directories(${PROJECT_SOURCE_DIR}/SuiteSparse/KLU/templates)
include_directories(${INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_SOURCE_DIR}/src)

######START:Make SuiteSParse#############
add_custom_command(
  COMMAND echo "Compiling SuiteSparse Library"
  COMMAND cd SuiteSparse && make -j5 library
  OUTPUT  ${PARENT_SOURCE_DIR}/SuiteSparse/KLU/lib/libklu.a
)

ADD_CUSTOM_TARGET(suitesparse ALL DEPENDS ${PARENT_SOURCE_DIR}/SuiteSparse/KLU/lib/libklu.a)
######END:Make SuiteSParse#############

#list all source files here
ADD_EXECUTABLE(blockklu main.cpp)
set(LIBS  ${LIBS} ${LIBSPARSE_LIBS}   ) 
target_link_libraries(blockklu ${LIBS})

#need to link to some other libraries ? just add them here
#TARGET_LINK_LIBRARIES(blockklu BlockKLu)
 
