 
#include <RcppParallel.h>
using namespace RcppParallel;
#include <RcppArmadillo.h>

#include "S6M_SellingReduce.h"
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Reducer.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
#include "D7Dist.h"
#include "CS6Dist.h"
#include "CS6Dist_func.h"
#include <cmath>
#include <stdlib.h>
#include <ctype.h>

// *****************************************************************************
G6 makeprimredcell( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma )
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 g6redprimcell;
    double g6primcell[6];
    double dg6redprimcell[6];
    Mat66 mc;
    Mat66 m;
    Mat66 dm;
    G6 primcell;
    LRL_Cell rawcell(a,b,c, alpha,beta,gamma);
    int ii;
    bool ret;
    int reduced;
    if (testlattice.size()< 1) {
        latsym = "P";
    } else {
        latsym = testlattice.substr(0,1);
    }
    clatsym = latsym[0];
    switch (clatsym) {
        case 'P':
        case 'p':
        case 'A':
        case 'a':
        case 'B':
        case 'b':
        case 'C':
        case 'c':
        case 'I':
        case 'i':
        case 'F':
        case 'f':
        case 'R':
        case 'r':
        case 'H':
        case 'h':
            CS6M_CelltoG6(rawcell,v6cell);
            CS6M_LatSymMat66(v6cell,clatsym,mc,primcell);
            break;
        case 'V':
        case 'v':
	    primcell[0] = a;
	    primcell[1] = b;
	    primcell[2] = c;
	    primcell[3] = alpha;
	    primcell[4] = beta;
	    primcell[5] = gamma;
            break;
        case 'D':
        case 'd':
           primcell[0] = a;
           primcell[1] = b;
           primcell[2] = c;
           primcell[3] = beta-b-c;
           primcell[4] = gamma-a-c;
           primcell[5] = alpha-c-primcell[3]-primcell[4] ;
           break;
        case 'S':
        case 's':
           primcell[3] = 2.*a;
           primcell[4] = 2.*b;
           primcell[5] = 2.*c;
           primcell[0] = -alpha-c-b;
           primcell[1] = -beta-c-a;
           primcell[2] = -gamma-b-a;
           break;
        default:
            /* Rprintf("Unrecognized lattice symbol %s treated as P\n",testlattice.c_str()); */
            latsym = "P";
            clatsym='P';
            CS6M_CelltoG6(rawcell,v6cell);
            CS6M_LatSymMat66(v6cell,clatsym,mc,primcell);
            break;
    }
    reduced=0;
    g6primcell[0]=primcell[0];   
    g6primcell[1]=primcell[1];   
    g6primcell[2]=primcell[2];   
    g6primcell[3]=primcell[3];   
    g6primcell[4]=primcell[4];   
    g6primcell[5]=primcell[5];
    CS6M_G6Reduce(g6primcell,dg6redprimcell,reduced);
    if (reduced) {
      g6redprimcell = G6(dg6redprimcell);
    } else {
      dg6redprimcell[0]=dg6redprimcell[1]=dg6redprimcell[2]
        =dg6redprimcell[3]=dg6redprimcell[4]=dg6redprimcell[5]=0;
      g6redprimcell = G6(dg6redprimcell);
    }
    return g6redprimcell;
}



extern "C" SEXP rcpp_cs6dist_in_g6 ( SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_, 
                              SEXP alpha1_, SEXP beta1_, SEXP gamma1_,
			      SEXP lat2_, SEXP a2_, SEXP b2_, SEXP c2_, 
                              SEXP alpha2_, SEXP beta2_, SEXP gamma2_) {
    G6 prim1, prim2;
    double dprim1[6];
    double dprim2[6];
    size_t ii;

    std::string lat1 = Rcpp::as<std::string>(lat1_);
    double a1        = Rcpp::as<double>(a1_);
    double b1        = Rcpp::as<double>(b1_);
    double c1        = Rcpp::as<double>(c1_);
    double alpha1    = Rcpp::as<double>(alpha1_);
    double beta1     = Rcpp::as<double>(beta1_);
    double gamma1    = Rcpp::as<double>(gamma1_);
    std::string lat2 = Rcpp::as<std::string>(lat2_);
    double a2        = Rcpp::as<double>(a2_);
    double b2        = Rcpp::as<double>(b2_);
    double c2        = Rcpp::as<double>(c2_);
    double alpha2    = Rcpp::as<double>(alpha2_);
    double beta2     = Rcpp::as<double>(beta2_);
    double gamma2    = Rcpp::as<double>(gamma2_);
    prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    prim2 = makeprimredcell(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    /* Rprintf("dprim1:  [%g, %g, %g, %g, %g, %g]\n", dprim1[0], dprim1[1], dprim1[2], dprim1[3], dprim1[4], dprim1[5]); */
    /* Rprintf("dprim2:  [%g, %g, %g, %g, %g, %g]\n", dprim2[0], dprim2[1], dprim2[2], dprim2[3], dprim2[4], dprim2[5]); */
    return Rcpp::wrap( 0.1*std::sqrt(CS6Dist_in_G6(dprim1,dprim2)));
}


