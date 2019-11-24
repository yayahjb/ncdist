######################################################################
#  CMakeLists.txt - cmake build file for make to create ncdist       #
#                                                                    #
# Version 1.1.0 07 Jul 2019                                          #
#                                                                    #
#         Herbert J. Bernstein (yayahjb@gmail.com)                   #
#         Lawrence C Andrews                                         #
#                                                                    #
# (C) Copyright 2016 - 2018 Herbert J. Bernstein, Lawrence C. Andrews#
#                                                                    #
######################################################################

######################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE NCDIST PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE NCDIST API UNDER THE TERMS  #
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
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

set(DEPENDENCIES
B4.cpp          B4.h
C3.cpp          C3.h
CellInputData.h CellInputData.cpp
D7.cpp          D7.h
D7.h            D7.cpp 
D7_Boundary.h   D7_Boundary.cpp
D7_BoundaryList.h D7_BoundaryList.cpp 
D7_Subboundary.h D7_Subboundary.cpp 
D7_ClassifySubboundaries.h D7_ClassifySubboundaries.cpp
G6.cpp  G6.h
inverse.cpp     inverse.h
LatticeConverter.h LatticeConverter.cpp
LRL_Cell.h      LRL_Cell.cpp  
LRL_Cell_Degrees.h   
LRL_Cell_Degrees.cpp  
LRL_CoordinateConversionMatrices.h LRL_CoordinateConversionMatrices.cpp  
LRL_StringTools.h  LRL_StringTools.cpp  
MatD7.h         MatD7.cpp  
MatG6.h         MatG6.cpp  
MatS6.h         MatS6.cpp  
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Mat66.cpp       Mat66.h
MaximaTools.h   MaximaTools.cpp
PrintTable.h    PrintTable.cpp
ProjectorTools.h ProjectorTools.cpp 
RandTools.h     RandTools.cpp 
ReadCellData.h  ReadCellData.cpp
Reducer.cpp     Reducer.h
S6.cpp          S6.h
S6Dist.cpp      S6Dist.h
Selling.cpp     Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp S6Dist_func.h
Vec_N_Tools.cpp Vec_N_Tools.h
VecN.cpp        VecN.h
VectorTools.cpp VectorTools.h
vector_3d.cpp   vector_3d.h)


set(CS6DIST_APP_SOURCES
cs6dist_app.cpp
CS6Dist_func.cpp	CS6Dist_func.h
CS6Dist.c CS6Dist.h
Delone.cpp              Delone.h
${DEPENDENCIES}
)

set(CS6DIST_DIST_SOURCES
cs6dist_dist.cpp
CS6Dist_func.cpp	CS6Dist_func.h
CS6Dist.c CS6Dist.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)


set(CS6DIST_MAT_SOURCES
cs6dist_mat.cpp
CS6Dist_func.cpp	CS6Dist_func.h
CS6Dist.c CS6Dist.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)

set(CS6_S6_TEST_SOURCES
cs6_s6_test.cpp
CS6Dist_func.cpp	CS6Dist_func.h
CS6Dist.c CS6Dist.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)


set(S6DIST_APP_SOURCES
s6dist_app.cpp
S6Dist_func.cpp	S6Dist_func.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)


set(MINIMAL_S6DIST_SOURCES
cqrlib.c        cqrlib.h
MatN.cpp        MatN.h
S6.cpp          S6.h
Selling.cpp     Selling.h
VecN.cpp        VecN.h
S6Dist.cpp      S6Dist.h
Vec_N_Tools.cpp Vec_N_Tools.h
inverse.cpp     inverse.h
MinimalS6Dist.cpp
Delone.cpp      Delone.h
S6.cpp          S6.h
triple.h
vector_3d.cpp   vector_3d.h
rhrand.h
${DEPENDENCIES}
)


set(NCDIST_SOURCES
NCDist.h    ncdist.cpp
Delone.h    Delone.cpp
${DEPENDENCIES}
)

set (RCPP_CS6DIST_SOURCES
rcpp_cs6dist.cpp
CS6Dist.h CS6Dist.c
S6Dist_func.cpp S6Dist_func.h
Delone.cpp Delone.h
${DEPENDENCIES}
)

set (RCPP_CS6DIST_IN_G6_SOURCES
rcpp_cs6dist_in_g6.cpp
CS6Dist.h CS6Dist.c
S6Dist_func.cpp S6Dist_func.h
Delone.cpp Delone.h
${DEPENDENCIES}
)

set(RCPP_D7DIST_SOURCES
rcpp_d7dist.cpp
D7Dist.c        D7Dist.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)


set(RCPP_NCDIST_SOURCES
rcpp_ncdist.cpp 
NCDist.h
Delone.cpp      Delone.h
${DEPENDENCIES}
)

set (RCPP_S6DIST_SOURCES
rcpp_s6dist.cpp
S6Dist_func.cpp  S6Dist_func.h
Delone.cpp Delone.h
${DEPENDENCIES}
)

set(D7DIST_SOURCES
d7dist.cpp
D7Dist.c        D7Dist.h
Delone.cpp Delone.h
${DEPENDENCIES}
)


set(FOLLOWER_SOURCES
D7Dist.c D7Dist.h
Delone.cpp Delone.h
CreateFileName.h CreateFileName.cpp
MapBoundaryStrings2Colors.h MapBoundaryStrings2Colors.cpp
FileOperations.cpp	FileOperations.h
FileWriter.h
Follow.h 
Follower.h 
FollowerConstants.cpp	FollowerConstants.h
FollowerIO.cpp	FollowerIO.h
FollowerTools.cpp	FollowerTools.h
Follower_main.cpp
LinearAxis.cpp LinearAxis.h
ReadGlobalData.cpp 
NCDist.c NCDist.h
${DEPENDENCIES}
)


add_library(R_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(R_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${R_LIBRARY} )
add_library(RBLAS_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RBLAS_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RBLAS_LIBRARY} )
add_library(RLAPACK_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RLAPACK_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RLAPACK_LIBRARY} )
add_library(rcpp_ncdist SHARED ${RCPP_NCDIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_d7dist SHARED ${RCPP_D7DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_cs6dist SHARED ${RCPP_CS6DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_cs6dist_in_g6 SHARED ${RCPP_CS6DIST_IN_G6_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_s6dist SHARED ${RCPP_S6DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
target_link_libraries(rcpp_ncdist)
add_executable(ncdist ${NCDIST_SOURCES} )
target_link_libraries(ncdist ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(rcpp_d7dist)
target_link_libraries(rcpp_cs6dist)
target_link_libraries(rcpp_cs6dist_in_g6)
add_executable(d7dist ${D7DIST_SOURCES} )
target_link_libraries(d7dist ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6dist_app ${CS6DIST_APP_SOURCES} )
target_link_libraries(cs6dist_app ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6dist_mat ${CS6DIST_MAT_SOURCES} )
target_link_libraries(cs6dist_mat ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6_s6_test ${CS6_S6_TEST_SOURCES} )
target_link_libraries(cs6_s6_test ${CMAKE_THREAD_LIBS_INIT} )
add_executable(s6dist_app ${S6DIST_APP_SOURCES} )
target_link_libraries(s6dist_app ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6dist_dist ${CS6DIST_DIST_SOURCES} )
target_link_libraries(cs6dist_dist ${CMAKE_THREAD_LIBS_INIT} )
add_executable(minimal_s6dist ${MINIMAL_S6DIST_SOURCES} )
target_link_libraries(minimal_s6dist ${CMAKE_THREAD_LIBS_INIT} )

add_executable(Follower ${FOLLOWER_SOURCES} )
target_link_libraries(Follower ${CMAKE_THREAD_LIBS_INIT} )


install(TARGETS cs6dist_app DESTINATION bin)
install(TARGETS cs6dist_mat DESTINATION bin)
install(TARGETS cs6_s6_test DESTINATION bin)
install(TARGETS cs6dist_dist DESTINATION bin)
install(TARGETS minimal_s6dist DESTINATION bin)
install(TARGETS ncdist DESTINATION bin)
install(TARGETS d7dist DESTINATION bin)
install(TARGETS Follower DESTINATION bin)
install(TARGETS s6dist_app DESTINATION bin)
install(TARGETS rcpp_cs6dist DESTINATION lib)
install(TARGETS rcpp_cs6dist_in_g6 DESTINATION lib)
install(TARGETS rcpp_d7dist DESTINATION lib)
install(TARGETS rcpp_ncdist DESTINATION lib)
install(TARGETS rcpp_s6dist DESTINATION lib)
install(FILES license.txt DESTINATION share/ncdist/)

