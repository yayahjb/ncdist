######################################################################
#  Makefile - makefile to create sfdist                              #
#                                                                    #
# Version 1.0.1 01 April 2018                                        #
#                                                                    #
#         Herbert J. Bernstein (yayahjb@gmail.com)                   #
#         Lawrence C Andrews                                         #
#                                                                    #
# (C) Copyright 2016 - 2018 Herbert J. Bernstein, Lawrence C. Andrews#
#                                                                    #
######################################################################

######################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE sfdist PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE sfdist API UNDER THE TERMS  #
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

ifeq ($(OS),Windows_NT)
    uname_S := Windows
else
    uname_S := $(shell uname -s)
endif

ifeq ($(uname_S), Darwin)
  RCPP_HEADERS ?= /Library/Frameworks/R.framework/Versions/3.3/Resources/library/Rcpp/include
  RCPPARMA_HEADERS ?=  /Library/Frameworks/R.framework/Versions/3.3/Resources/library/RcppArmadillo/include
  RCPPPARA_HEADERS ?= /Library/Frameworks/R.framework/Versions/3.3/Resources/library/RcppParallel/include
  RPATH_HEADERS ?= /Library/Frameworks/R.framework/Versions/3.3/Headers
  RPATH_LIBRARIES ?= /Library/Frameworks/R.framework/Versions/3.3/Resources/lib
else
  RCPP_HEADERS ?= /usr/local/lib/R/site-library/Rcpp/include
  RCPPARMA_HEADERS ?= /usr/local/lib/R/site-library/RcppArmadillo/include
  RCPPPARA_HEADERS ?= /usr/local/lib/R/site-library/RcppParallel/include
  RPATH_HEADERS ?= /usr/share/R/include
  RPATH_LIBRARIES ?= /usr/lib/R/lib
endif

CFLAGS ?= -g -O3 -fPIC -fopenmp
CXXFLAGS ?= -g -O3 -std=gnu++11 -fPIC -fopenmp
#CFLAGS ?= -g -O0 -fPIC -fopenmp
#CXXFLAGS ?= -g -O0 -std=gnu++11 -fPIC -fopenmp


LIBSOURCES = \
	B4.cpp \
	C3.cpp \
	CellInputData.cpp \
	D7.cpp \
	D7_Boundary.cpp \
	D7_BoundaryList.cpp \
	D7_Subboundary.cpp \
	D7_ClassifySubboundaries.cpp\
	G6.cpp \
	inverse.cpp \
	LatticeConverter.cpp \
	LRL_Cell.cpp \
	LRL_Cell_Degrees.cpp \
	LRL_CoordinateConversionMatrices.cpp \
	LRL_StringTools.cpp \
	MatD7.cpp \
	MatG6.cpp \
	MatS6.cpp \
	MatN.cpp \
	MatMN.cpp \
	Mat66.cpp \
	MaximaTools.cpp \
	PrintTable.cpp \
	ProjectorTools.cpp \
	RandTools.cpp \
	ReadCellData.cpp \
	Reducer.cpp \
	S6.cpp \
	S6Dist.cpp \
	Selling.cpp \
	S6M_SellingReduce.h \
	VectorTools.h \
	VecN.h VecN.cpp \
	Vec_N_Tools.cpp \
	VectorTools.cpp \
	vector_3d.cpp 


DEPENDENCIES =  \
	B4.h B4.cpp \
	CellInputData.h CellInputData.cpp \
	D7.h D7.cpp \
	D7_Boundary.h D7_Boundary.cpp \
	D7_BoundaryList.h D7_BoundaryList.cpp \
	D7_Subboundary.h D7_Subboundary.cpp \
	D7_ClassifySubboundaries.h D7_ClassifySubboundaries.cpp\
	G6.h G6.cpp \
	inverse.h inverse.cpp \
	LatticeConverter.h LatticeConverter.cpp \
	LRL_Cell.h LRL_Cell.cpp \
	LRL_Cell_Degrees.h  \
	LRL_Cell_Degrees.cpp \
	LRL_CoordinateConversionMatrices.h LRL_CoordinateConversionMatrices.cpp \
	LRL_StringTools.h  LRL_StringTools.cpp \
	MatD7.h MatD7.cpp \
	MatG6.h MatG6.cpp \
	MatS6.h MatS6.cpp \
	MatN.h MatN.cpp \
	MatMN.h MatMN.cpp \
	Mat66.h Mat66.cpp \
	MaximaTools.h MaximaTools.cpp \
	PrintTable.h PrintTable.cpp \
	ProjectorTools.h ProjectorTools.cpp \
	RandTools.h RandTools.cpp \
	ReadCellData.h ReadCellData.cpp \
	Reducer.h Reducer.cpp \
	S6.h S6.cpp \
	S6M_SellingReduce.h \
	Vec_N_Tools.h VectorTools.h \
	VecN.h VecN.cpp \
	Vec_N_Tools.cpp \
	VectorTools.cpp \
	vector_3d.cpp 



all:  ncdist ncdist_mat D7Test Follower \
	rcpp_ncdist.so rcpp_d7dist.so rcpp_s6dist.so rcpp_cs6dist.so \
	cs6dist_app  s6dist_app cs6dist_mat cs6dist_dist cs6_s6_test

CS6Dist_func_.o:  CS6Dist_func.cpp S6.h C3.h S6Dist.h Selling.h \
	S6Dist_func.h Reducer.h Delone.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h \
	CS6Dist.h VecN.h Vec_N_Tools.h 
	g++ $(CXXFLAGS) -c -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	CS6Dist_func.cpp -o CS6Dist_func_.o

cs6dist_app:  CS6Dist_func_.o cs6dist_app.cpp \
	$(DEPENDENCIES) Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o cs6dist_app  cs6dist_app.cpp CS6Dist_func_.o \
	-I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES)  Delone.cpp \
	-lpthread


cs6_s6_test:  CS6Dist_func_.o S6Dist_func_.o cs6_s6_test.cpp \
	$(DEPENDENCIES) Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o cs6_s6_test  cs6_s6_test.cpp CS6Dist_func_.o S6Dist_func_.o \
	-I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES) Delone.cpp \
	-lpthread


cs6dist_mat:  CS6Dist_func_.o cs6dist_mat.cpp \
	$(DEPENDENCIES)  Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o cs6dist_mat  cs6dist_mat.cpp CS6Dist_func_.o \
	-I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES)  Delone.cpp \
	-lpthread

cs6dist_dist:  CS6Dist_func_.o cs6dist_dist.cpp \
	$(DEPENDENCIES)  Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o cs6dist_dist  cs6dist_dist.cpp CS6Dist_func_.o \
	-I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES)  Delone.cpp \
	-lpthread

S6Dist_func_.o:  S6Dist_func.cpp S6.h C3.h S6Dist.h Selling.h \
	S6Dist_func.h Reducer.h Delone.h LRL_Cell_Degrees.h NCDist.h VecN.h Vec_N_Tools.h 
	g++ $(CXXFLAGS) -c -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) S6Dist_func.cpp -o S6Dist_func_.o

s6dist_app:  S6Dist_func_.o s6dist_app.cpp \
	$(DEPENDENCIES)  Delone.h Delone.cpp 
	g++ $(CXXFLAGS) -o s6dist_app  s6dist_app.cpp S6Dist_func_.o \
	-I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES)  Delone.cpp \
	-lpthread

ncdist_mat_.o:  ncdist_mat.cpp Reducer.h Delone.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h VecN.h Vec_N_Tools.h
	g++ $(CXXFLAGS) -c -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) ncdist_mat.cpp -o ncdist_mat_.o

ncdist_mat:  ncdist_mat_.o \
	$(DEPENDENCIES) Delone.h Delone.cpp 
	g++ $(CXXFLAGS) -o ncdist_mat ncdist_mat_.o -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES) Delone.cpp \
	-lpthread

ncdist_.o:  ncdist.cpp Reducer.h Delone.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h VecN.h Vec_N_Tools.h
	g++ $(CXXFLAGS) -c -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) ncdist.cpp -o ncdist_.o

ncdist:  ncdist_.o \
	$(DEPENDENCIES) Delone.h Delone.cpp  
	g++ $(CXXFLAGS) -o ncdist ncdist_.o -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES)Delone.cpp \
	-lpthread

d7dist_.o:  d7dist.cpp Reducer.h Delone.h LRL_Cell.h  LRL_Cell_Degrees.h NCDist.h VecN.h Vec_N_Tools.h
	g++ $(CXXFLAGS) -c -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) d7dist.cpp -o d7dist_.o

d7dist:  d7dist_.o \
	$(DEPENDENCIES) Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o d7dist d7dist_.o -I $(RCPPARMA_HEADERS) -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES) \
	-lpthread

rcpp_cs6dist_.o:  rcpp_cs6dist.cpp G6.h D7.h S6.h C3.h Delone.h Reducer.h LRL_Cell.h LRL_Cell_Degrees.h \
	D7Dist.h S6Dist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c rcpp_cs6dist.cpp -o rcpp_cs6dist_.o

rcpp_s6dist_.o:  rcpp_s6dist.cpp G6.h D7.h S6.h C3.h Delone.h Reducer.h LRL_Cell_Degrees.h \
	D7Dist.h S6Dist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c rcpp_s6dist.cpp -o rcpp_s6dist_.o

rcpp_ncdist_.o:  rcpp_ncdist.cpp Reducer.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c rcpp_ncdist.cpp -o rcpp_ncdist_.o

rcpp_d7dist_.o:  rcpp_d7dist.cpp Reducer.h LRL_Cell.h LRL_Cell_Degrees.h D7Dist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c rcpp_d7dist.cpp -o rcpp_d7dist_.o

Reducer.o:  Reducer.cpp Reducer.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c Reducer.cpp -o Reducer.o


Delone.o:  Delone.cpp Delone.h  Reducer.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c Delone.cpp  -o Delone.o

LRL_Cell.o:  LRL_Cell.cpp LRL_Cell_Degrees.cpp Reducer.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c LRL_Cell.cpp -o LRL_Cell.o

LRL_Cell_Degrees.o:  LRL_Cell.cpp LRL_Cell_Degrees.cpp Reducer.h LRL_Cell.h LRL_Cell_Degrees.h NCDist.h
	g++ $(CXXFLAGS)  -I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP_HEADERS) \
		-I$(RCPPPARA_HEADERS) \
		-I$(RCPPARMA_HEADERS) -c LRL_Cell_Degrees.cpp -o LRL_Cell_Degrees.o

rcpp_s6dist.so:	rcpp_s6dist_.o Reducer.o LRL_Cell.o \
	$(DEPENDENCIES) \
	S6Dist_func.cpp
	g++ $(CXXFLAGS) -shared -o rcpp_s6dist.so -I $(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC \
        -I $(RCPP_HEADERS) -I $(RCPPPARA_HEADERS) -I$(RCPPARMA_HEADERS) rcpp_s6dist.cpp  \
	$(LIBSOURCES) \
	S6Dist_func.cpp -L$(RPATH_LIBRARIES) -lR -lblas -llapack -lpthread

rcpp_cs6dist.so: rcpp_cs6dist_.o Reducer.o LRL_Cell.o LRL_Cell_Degrees.o \
	CS6Dist.h $(DEPENDENCIES) 
	g++ $(CXXFLAGS) -shared -o rcpp_cs6dist.so rcpp_cs6dist_.o -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES) \
	-L$(RPATH_LIBRARIES) -lR -lblas -llapack -lpthread

rcpp_ncdist.so:	rcpp_ncdist_.o Reducer.o LRL_Cell.o LRL_Cell_Degrees.o
	g++ $(CXXFLAGS) -shared -o rcpp_ncdist.so rcpp_ncdist_.o -I $(RCPPPARA_HEADERS) \
	$(LIBSOURCES) \
	-L$(RPATH_LIBRARIES) -lR -lblas -llapack -lpthread

rcpp_d7dist.so:	rcpp_d7dist_.o Reducer.o LRL_Cell.o LRL_Cell_Degrees.o \
	$(DEPENDENCIES)
	g++ $(CXXFLAGS) -shared -o rcpp_d7dist.so \
	-I$(RPATH_HEADERS) -DNDEBUG  -fpic  -O2 -fPIC -I$(RCPP_HEADERS) \
	-I$(RCPPPARA_HEADERS)  -I$(RCPPARMA_HEADERS) rcpp_d7dist.cpp \
	$(LIBSOURCES) \
	-L$(RPATH_LIBRARIES) -lR -lblas -llapack -lpthread

NCDist.o:  NCDist.c NCDist.h
	gcc $(CFLAGS) -c NCDist.c

D7Dist.o:  D7Dist.c D7Dist.h
	gcc $(CFLAGS) -c D7Dist.c

D7Test: D7Test.cpp D7Dist.h NCDist.o  \
	$(DEPENDENCIES) Delone.h Delone.cpp
	g++ $(CXXFLAGS) -o D7Test  D7Test.cpp  NCDist.o Delone.cpp \
	$(LIBSOURCES) \
	-lpthread

Follower: \
	$(DEPENDENCIES) \
	Delone.cpp \
	CreateFileName.h CreateFileName.cpp \
	MapBoundaryStrings2Colors.h MapBoundaryStrings2Colors.cpp \
	FileOperations.cpp FileOperations.cpp \
	FollowerConstants.h FollowerConstants.cpp \
	FollowerIO.h FollowerIO.cpp \
	FollowerTools.h FollowerTools.cpp \
	Follower_main.cpp \
	NCDist.o \
	FileOperations.h FileWriter.h Follow.h Follower.h \
	FollowerConstants.h FollowerIO.h FollowerTools.h \
	LinearAxis.h LinearAxis.cpp \
	ReadGlobalData.h ReadGlobalData.cpp
	g++ $(CXXFLAGS) -o Follower \
	$(LIBSOURCES) \
	Delone.cpp \
	CreateFileName.cpp \
	MapBoundaryStrings2Colors.cpp \
	FileOperations.cpp \
	FollowerConstants.cpp FollowerIO.cpp FollowerTools.cpp \
	LinearAxis.cpp \
	Follower_main.cpp \
	ReadGlobalData.cpp \
	NCDist.o 

MinimalS6Dist:  cqrlib.c     MatN.cpp           S6.cpp C3.cpp      Selling.cpp  \
	VecN.h	cqrlib.h     MatN.h             S6Dist.cpp  Selling.h    \
	Vec_N_Tools.cpp	     inverse.cpp        MinimalS6Dist.cpp \
	S6Dist.h             TNear.h            Vec_N_Tools.h \
	inverse.h            S6.h C3.h               triple.h     \
	vector_3d.cpp        rhrand.h           VecN.cpp    vector_3d.h
	g++ $(CXXFLAGS) -o   MinimalS6Dist \
	MinimalS6Dist.cpp \
	inverse.cpp          S6Dist.cpp         VecN.cpp    vector_3d.cpp \
	MatN.cpp             S6.cpp C3.cpp             Selling.cpp Vec_N_Tools.cpp 

clean:
	-@rm -rf *.o
	-@rm -rf *.so
	-@rm -rf D7Test
	-@rm -rf Follower

distclean:  clean
	-@rm -rf bin
	-@rm -rf build
	-@rm -rf CS6Dist_func_.o
	-@rm -rf cs6dist_app
	-@rm -rf cs6_s6_test
	-@rm -rf cs6dist_mat
	-@rm -rf cs6dist_dist
	-@rm -rf d7dist
	-@rm -rf d7dist_.o
	-@rm -rf Delone.o
	-@rm -rf D7Test
	-@rm -rf Follower
	-@rm -rf lib
	-@rm -rf LRL_Cell.o
	-@rm -rf LRL_Cell_Degrees.o
	-@rm -rf MinimalS6Dist
	-@rm -rf ncdist
	-@rm -rf ncdist_mat_.o
	-@rm -rf ncdist_mat
	-@rm -rf ncdist_.o
	-@rm -rf ncdist
	-@rm -rf rcpp_cs6dist_.o
	-@rm -rf rcpp_s6dist_.o
	-@rm -rf rcpp_ncdist_.o
	-@rm -rf rcpp_d7dist_.o
	-@rm -rf Reducer.o
	-@rm -rf S6Dist_func_.o
	-@rm -rf s6dist_app
	-@rm -rf rcpp_s6dist.so
	-@rm -rf rcpp_cs6dist.so
	-@rm -rf rcpp_ncdist.so
	-@rm -rf rcpp_d7dist.so
