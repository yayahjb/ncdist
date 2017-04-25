RCPP ?= /usr/lib/R/site-library/Rcpp
RCPPARMA ?= /usr/lib/R/site-library/RcppArmadillo
RCPPPARA ?= /usr/local/lib/R/site-library/RcppParallel/
RPATH ?= /usr/share/R

ncdist.o:  ncdist.cpp Reducer.h Cell.h NCDist.h
	g++ -g -O3 -fopenmp -c -I $(RCPPARMA)/include -I $(RCPPPARA)/include ncdist.cpp

ncdist:  ncdist.o Reducer.cpp Cell.cpp BasicDistance.cpp BasicDistance.h
	g++ -g -O3 -fopenmp -o ncdist ncdist.o -I $(RCPPARMA)/include -I $(RCPPPARA)/include Reducer.cpp Cell.cpp BasicDistance.cpp -lpthread

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


D7Test:  Cell.cpp CellInputData.cpp D6.cpp D7.cpp D7Test.cpp Delaunay.cpp DelaunayTetrahedron.cpp \
	G6.cpp LatticeConverter.cpp Mat66.cpp Mat7.cpp MatMN.cpp MatN.cpp MaximaTools.cpp PrintTable.cpp \
	ProjectorTools.cpp RandTools.cpp ReadCellData.cpp Reducer.cpp StringTools.cpp VecN.cpp \
	Vec_N_Tools.cpp VectorTools.cpp inverse.cpp vector_3d.cpp D7dist.h NCDist.h
	g++ -g -O3 -fopenmp -o D7Test  \
	Cell.cpp CellInputData.cpp D6.cpp D7.cpp D7Test.cpp Delaunay.cpp DelaunayTetrahedron.cpp \
	G6.cpp LatticeConverter.cpp Mat66.cpp Mat7.cpp MatMN.cpp MatN.cpp MaximaTools.cpp PrintTable.cpp \
	ProjectorTools.cpp RandTools.cpp ReadCellData.cpp Reducer.cpp StringTools.cpp VecN.cpp \
	Vec_N_Tools.cpp VectorTools.cpp inverse.cpp vector_3d.cpp -lpthread

