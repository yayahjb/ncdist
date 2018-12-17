######################################################################
#  CMakeLists.txt - cmake build file for make to create ncdist       #
#                                                                    #
# Version 1.0.3 16 Dec 2018                                          #
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

set(CS6DIST_APP_SOURCES
cs6dist_app.cpp
CS6Dist_func.cpp	CS6Dist_func.h
Cell.cpp        Cell.h
CS6Dist.h
D6.cpp	D6.h
D7.cpp	D7.h
Delone.cpp	Delone.h
DeloneTetrahedron.cpp	DeloneTetrahedron.h
G6.cpp	G6.h
inverse.cpp 	inverse.h
MatN.cpp	MatN.h
MatMN.cpp	MatMN.h
Mat66.cpp	Mat66.h
NCDist.h
Reducer.cpp	Reducer.h
S6.cpp          S6.h
S6Dist.cpp	S6Dist.h
Selling.cpp	Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp	S6Dist_func.h
Vec_N_Tools.cpp	Vec_N_Tools.h
VecN.cpp	VecN.h
VectorTools.cpp	VectorTools.h
vector_3d.cpp	vector_3d.h)

set(CS6DIST_DIST_SOURCES
cs6dist_dist.cpp
CS6Dist_func.cpp	CS6Dist_func.h
Cell.cpp        Cell.h
CS6Dist.h
D6.cpp	D6.h
D7.cpp	D7.h
Delone.cpp	Delone.h
DeloneTetrahedron.cpp	DeloneTetrahedron.h
G6.cpp	G6.h
inverse.cpp 	inverse.h
MatN.cpp	MatN.h
MatMN.cpp	MatMN.h
Mat66.cpp	Mat66.h
NCDist.h
Reducer.cpp	Reducer.h
S6.cpp          S6.h
S6Dist.cpp	S6Dist.h
Selling.cpp	Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp	S6Dist_func.h
Vec_N_Tools.cpp	Vec_N_Tools.h
VecN.cpp	VecN.h
VectorTools.cpp	VectorTools.h
vector_3d.cpp	vector_3d.h)


set(CS6DIST_MAT_SOURCES
cs6dist_mat.cpp
CS6Dist_func.cpp	CS6Dist_func.h
Cell.cpp        Cell.h
CS6Dist.h
D6.cpp	D6.h
D7.cpp	D7.h
Delone.cpp	Delone.h
DeloneTetrahedron.cpp	DeloneTetrahedron.h
G6.cpp	G6.h
inverse.cpp 	inverse.h
MatN.cpp	MatN.h
MatMN.cpp	MatMN.h
Mat66.cpp	Mat66.h
NCDist.h
Reducer.cpp	Reducer.h
S6.cpp          S6.h
S6Dist.cpp	S6Dist.h
Selling.cpp	Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp	S6Dist_func.h
Vec_N_Tools.cpp	Vec_N_Tools.h
VecN.cpp	VecN.h
VectorTools.cpp	VectorTools.h
vector_3d.cpp	vector_3d.h)


set(S6DIST_APP_SOURCES
s6dist_app.cpp
S6Dist_func.cpp	S6Dist_func.h
Cell.cpp        Cell.h
CS6Dist.h
D6.cpp	D6.h
D7.cpp	D7.h
Delone.cpp	Delone.h
DeloneTetrahedron.cpp	DeloneTetrahedron.h
G6.cpp	G6.h
inverse.cpp 	inverse.h
MatN.cpp	MatN.h
MatMN.cpp	MatMN.h
Mat66.cpp	Mat66.h
NCDist.h
Reducer.cpp	Reducer.h
S6.cpp          S6.h
S6Dist.cpp	S6Dist.h
Selling.cpp	Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp	S6Dist_func.h
Vec_N_Tools.cpp	Vec_N_Tools.h
VecN.cpp	VecN.h
VectorTools.cpp	VectorTools.h
vector_3d.cpp	vector_3d.h)


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
S6.cpp          S6.h
triple.h
vector_3d.cpp   vector_3d.h
rhrand.h
)


set(NCDIST_SOURCES
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
NCDist.h    ncdist.cpp)

set (RCPP_CS6DIST_SOURCES
rcpp_cs6dist.cpp
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
S6.cpp          S6.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
D7Dist.h )

set(RCPP_D7DIST_SOURCES
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
S6.cpp          S6.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
D7Dist.h        rcpp_d7dist.cpp)


set(RCPP_NCDIST_SOURCES
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
S6.cpp          S6.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
NCDist.h        rcpp_ncdist.cpp)

set (RCPP_S6DIST_SOURCES
rcpp_s6dist.cpp
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
S6.cpp          S6.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
D7Dist.h )

set(D7DIST_SOURCES
Cell.cpp        Cell.h
D6.cpp          D6.h
D7.cpp          D7.h
Delone.cpp      Delone.h
DeloneTetrahedron.cpp
DeloneTetrahedron.h
G6.cpp          G6.h
inverse.cpp inverse.h
Mat66.cpp       Mat66.h
MatN.cpp        MatN.h
MatMN.cpp       MatMN.h
Reducer.cpp     Reducer.h
VecN.cpp        VecN.h
Vec_N_Tools.cpp Vec_N_Tools.h
vector_3d.cpp   vector_3d.h
VectorTools.cpp VectorTools.h
D7Dist.h    d7dist.cpp)


set(FOLLOWER_SOURCES
Cell.cpp        Cell.h
CellInputData.cpp	CellInputData.h 
cqrlib.c	cqrlib.h
CreateFileName.cpp	CreateFileName.h
CS6Dist.h
D6.cpp	D6.h
D7.cpp	D7.h
D7Dist.c	D7Dist.h	D7Dist_.hpp
Delone.cpp	Delone.h
DeloneTetrahedron.cpp	DeloneTetrahedron.h
FileOperations.cpp	FileOperations.h
FileWriter.h
Follow.h 
Follower.h 
FollowerConstants.cpp	FollowerConstants.h
FollowerIO.cpp	FollowerIO.h
FollowerTools.cpp	FollowerTools.h
Follower_main.cpp 
G6.cpp	G6.h
Glitch.h GlitchElement.h
inverse.cpp 	inverse.h
LinearAxis.cpp
MapBoundaryStrings2Colors.cpp	MapBoundaryStrings2Colors.h
Mat7.cpp	Mat7.h
MaximaTools.cpp	MaximaTools.h
MatN.cpp	MatN.h
MatMN.cpp	MatMN.h
Mat66.cpp	Mat66.h
NCDist.c	NCDist.h	NCDist_.hpp
NormedInt.h
ProjectorTools.cpp ProjectorTools.h
RandTools.cpp	RandTools.h
ReadCellData.cpp	ReadCellData.h 
ReadGlobalData.cpp	ReadGlobalData.h
Reducer.cpp	Reducer.h
ReportTools.cpp	ReportTools.h
Rn_Path.h 
rhrand.h
S6.cpp          S6.h
S6Dist.cpp	S6Dist.h
Selling.cpp	Selling.h
S6M_SellingReduce.h
S6Dist_func.cpp	S6Dist_func.h
SVG_CirclePlot.h	SVG_DistancePlot.h	SVG_Writer.h
StoreResults.h 
SuperFastHash.cpp	SuperFastHash.h 
TNear.h Theta.h
ToString.h
triple.h
Vec_N_Tools.cpp	Vec_N_Tools.h
VecN.cpp	VecN.h
VectorTools.cpp	VectorTools.h
vector_3d.cpp	vector_3d.h)


add_library(R_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(R_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${R_LIBRARY} )
add_library(RBLAS_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RBLAS_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RBLAS_LIBRARY} )
add_library(RLAPACK_shared_lib INTERFACE IMPORTED GLOBAL)
set_target_properties(RLAPACK_shared_lib PROPERTIES INTERFACE_LINK_LIBRARIES ${RLAPACK_LIBRARY} )
add_library(rcpp_ncdist SHARED ${RCPP_NCDIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_d7dist SHARED ${RCPP_D7DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_cs6dist SHARED ${RCPP_CS6DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
add_library(rcpp_s6dist SHARED ${RCPP_S6DIST_SOURCES} ${R_LIBRARY} ${RCPP_LIBRARY})
target_link_libraries(rcpp_ncdist)
add_executable(ncdist ${NCDIST_SOURCES} )
target_link_libraries(ncdist ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries(rcpp_d7dist)
target_link_libraries(rcpp_cs6dist)
add_executable(d7dist ${D7DIST_SOURCES} )
target_link_libraries(d7dist ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6dist_app ${CS6DIST_APP_SOURCES} )
target_link_libraries(cs6dist_app ${CMAKE_THREAD_LIBS_INIT} )
add_executable(cs6dist_mat ${CS6DIST_MAT_SOURCES} )
target_link_libraries(cs6dist_mat ${CMAKE_THREAD_LIBS_INIT} )
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
install(TARGETS cs6dist_dist DESTINATION bin)
install(TARGETS minimal_s6dist DESTINATION bin)
install(TARGETS ncdist DESTINATION bin)
install(TARGETS d7dist DESTINATION bin)
install(TARGETS Follower DESTINATION bin)
install(TARGETS s6dist_app DESTINATION bin)
install(TARGETS rcpp_cs6dist DESTINATION lib)
install(TARGETS rcpp_d7dist DESTINATION lib)
install(TARGETS rcpp_ncdist DESTINATION lib)
install(TARGETS rcpp_s6dist DESTINATION lib)
install(FILES license.txt DESTINATION share/ncdist/)

