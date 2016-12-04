RCPP ?= /usr/lib/R/site-library/Rcpp
RCPPARMA ?= /usr/lib/R/site-library/RcppArmadillo
RCPPPARA ?= /usr/local/lib/R/site-library/RcppParallel/
RPATH ?= /usr/share/R

ncdist.o:  ncdist.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp -c -I $(HOME)/usr/include ncdist.cpp

ncdist:  ncdist.o Reducer.cpp Cell.cpp BasicDistance.cpp BasicDistance.h
	g++ -g -O3 -fopenmp -o ncdist ncdist.o Reducer.cpp Cell.cpp BasicDistance.cpp -lpthread

rcpp_ncdist.o:  rcpp_ncdist.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP)/include \
		-I$(RCPPPARA)include \
		-I$(RCPPARMA)/include -c rcpp_ncdist.cpp -o rcpp_ncdist.o

Reducer.o:  Reducer.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP)include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c Reducer.cpp -o Reducer.o

Cell.o:  Cell.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP)include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c Cell.cpp -o Cell.o

BasicDistance.o:  BasicDistance.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp  -I$(RPATH)/include -DNDEBUG  -fpic  -O2 -fPIC \
	 -I$(RCPP)/include \
		-I$(RCPPPARA)/include \
		-I$(RCPPARMA)/include -c BasicDistance.cpp -o BasicDistance.o

rcpp_ncdist.so:	rcpp_ncdist.o Reducer.o Cell.o BasicDistance.o
	g++ -shared -o rcpp_ncdist.so rcpp_ncdist.o Reducer.o Cell.o BasicDistance.o -L/usr/lib/R/lib -lR -lblas -llapack -lpthread

