######################################################################
#  CMakeLists.txt - cmake build file for make to create ncdist       #
#                                                                    #
# Version 1.0.0 03 Dec 2016                                          #
#                                                                    #
#         Herbert J. Bernstein (yayahjb@gmail.com)                   #
#         Lawrence C Andrews                                         #
#                                                                    #
# (C) Copyright 2016 Herbert J. Bernstein, Lawrence C. Andrews       #
#                                                                    #
######################################################################

######################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE CBFLIB PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE CBFLIB API UNDER THE TERMS  #
# OF THE LGPL                                                        #
#                                                                    #
######################################################################

########################### GPL NOTICES ##############################
#                                                                    #
# This program is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU General Public License as     #
# published by the Free Software Foundation; either version 2 of     #
# (the License, or (at your option) any later version.               #
#                                                                    #
# This program is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      #
# GNU General Public License for more details.                       #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with this program; if not, write to the Free Software        #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           #
# 02111-1307  USA                                                    #
#                                                                    #
######################################################################

######################### LGPL NOTICES ###############################
#                                                                    #
# This library is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU Lesser General Public         #
# License as published by the Free Software Foundation; either       #
# version 2.1 of the License, or (at your option) any later version. #
#                                                                    #
# This library is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# Lesser General Public License for more details.                    #
#                                                                    #
# You should have received a copy of the GNU Lesser General Public   #
# License along with this library; if not, write to the Free         #
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    #
# MA  02110-1301  USA                                                #
#                                                                    #
######################################################################

cmake_minimum_required(VERSION 2.8)
project(ncdist C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(R REQUIRED)
find_package(Rcpp REQUIRED)
find_package(RcppArmadillo REQUIRED)
find_package(RcppParallel)

find_package(OpenMP)
if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")


find_package( Threads )

include_directories("${CMAKE_SOURCE_DIR}" "${R_INCLUDE_DIRS}"
"${RCPP_INCLUDE_DIRS}" "${RCPPARMA_INCLUDE_DIRS}" "${RCPPPARA_INCLUDE_DIRS}")

message("include RCPPPARA: "  ${RCPPPARA_INCLUDE_DIRS})

set(NCDIST_SOURCES
BasicDistance.cpp   BasicDistance.h     Cell.cpp    Cell.h
Reducer.cpp         Reducer.h           NCDist.h    ncdist.cpp)

set(RCPP_NCDIST_SOURCES
BasicDistance.cpp   BasicDistance.h     Cell.cpp    Cell.h
Reducer.cpp         Reducer.h           NCDist.h    rcpp_ncdist.cpp)


add_library(R_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(R_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${R_LIBRARY} )
add_library(RBLAS_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RBLAS_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RBLAS_LIBRARY} )
add_library(RLAPACK_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RLAPACK_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RLAPACK_LIBRARY} )
add_library(rcpp_ncdist SHARED ${RCPP_NCDIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
target_link_libraries(rcpp_ncdist R_shared_lib RBLAS_shared_lib RLAPACK_shared_lib)
add_executable(ncdist ${NCDIST_SOURCES} )
target_link_libraries(ncdist ${CMAKE_THREAD_LIBS_INIT} )

install(TARGETS ncdist DESTINATION bin)
install(TARGETS rcpp_ncdist DESTINATION lib)
install(FILES license.txt DESTINATION share/ncdist/)

