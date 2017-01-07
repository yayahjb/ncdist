######################################################################
#  Makefile - makefile to create sfdist                              #
#                                                                    #
# Version 1.0.0 07 January 2016                                      #
#                                                                    #
#         Herbert J. Bernstein (yayahjb@gmail.com)                   #
#         Lawrence C Andrews                                         #
#                                                                    #
# (C) Copyright 2016 - 2017 Herbert J. Bernstein, Lawrence C. Andrews#
#                                                                    #
######################################################################

######################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE sfdist PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE sfdist API UNDER THE TERMS  #
# OF THE sfdist                                                      #
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


PREFIX ?= $(CCP4)

RCPP ?= /usr/lib/R/site-library/Rcpp
RCPPARMA ?= /usr/lib/R/site-library/RcppArmadillo
RCPPPARA ?= /usr/local/lib/R/site-library/RcppParallel/
RPATH ?= /usr/share/R

bin:
	mkdir -p bin

lib:
	mkdir -p lib

build:
	mkdir -p build


all:	lib/librcpp_ncdist.so bin/ncdist

install:  all
	cp lib/librcpp_ncdist.so $(PREFIX)/lib/librcpp_ncdist.so
	chmod 755 $(PREFIX)/lib/librcpp_ncdist.so 
	cp bin/ncdist $(PREFIX)/bin/ncdist
	chmod 755 $(PREFIX)/bin/ncdist

clean:
	-@rm -f *.o

distclean:  clean
	-@rm -rf bin
	-@rm -rf lib
	-@rm -rf build

ncdist.o:  ncdist.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp -c -I $(RCPPARMA)/include -I $(RCPPPARA)/include ncdist.cpp

bin/ncdist:  bin ncdist.o Reducer.cpp Cell.cpp BasicDistance.cpp BasicDistance.h
	g++ -g -O3 -fopenmp -o bin/ncdist ncdist.o -I $(RCPPARMA)/include -I $(RCPPPARA)/include Reducer.cpp Cell.cpp BasicDistance.cpp -lpthread

rcpp_ncdist.o:  rcpp_ncdist.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -fPIC \
	 -I$(RCPP)/include \
		-I$(RCPPPARA)include \
		-I$(RCPPARMA)/include -c rcpp_ncdist.cpp -o rcpp_ncdist.o

Reducer.o:  Reducer.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -fPIC \
	 -I$(RCPP)include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c Reducer.cpp -o Reducer.o

Cell.o:  Cell.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic -fPIC \
	 -I$(RCPP)include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c Cell.cpp -o Cell.o

BasicDistance.o:  BasicDistance.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -fPIC \
	 	-I$(RCPP)/include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c BasicDistance.cpp -o BasicDistance.o

lib/librcpp_ncdist.so:	lib rcpp_ncdist.o Reducer.o Cell.o BasicDistance.o
	g++ -shared -o lib/librcpp_ncdist.so rcpp_ncdist.o Reducer.o Cell.o \
		BasicDistance.o -L/usr/lib/R/lib -lR -lblas -llapack -lpthread

