
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
#include <cmath>
#include <stdlib.h>
#include <ctype.h>


//*****************************************************************************
D7 makeprimredcell( std::string testlattice,
        double a, double b, double c, double alpha, double beta, double gamma )
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    D7 d7redprimcell;
    double dredprimcell[7];
    double g6primcell[6];
    double d7primcell[7];
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
    CS6M_G6toD7(g6primcell,d7primcell);
    CS6M_D7Reduce(d7primcell,dredprimcell,reduced);
    if (reduced) {
      d7redprimcell = D7(dredprimcell);
    } else {
      dredprimcell[0]=dredprimcell[1]=dredprimcell[2]=dredprimcell[3]
        =dredprimcell[4]=dredprimcell[5]=dredprimcell[6]=0.;
      d7redprimcell = D7(dredprimcell);
    }
    return d7redprimcell;
}



extern "C" SEXP rcpp_d7dist ( SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_, 
                              SEXP alpha1_, SEXP beta1_, SEXP gamma1_,
			      SEXP lat2_, SEXP a2_, SEXP b2_, SEXP c2_, 
                              SEXP alpha2_, SEXP beta2_, SEXP gamma2_) {
    G6 prim1, prim2;
    D7 d7prim1;
    D7 d7prim2;
    double dprim1[7];
    double dprim2[7];
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
    d7prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    d7prim2 = makeprimredcell(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    for (ii=0; ii < 7; ii++ ) dprim1[ii] = d7prim1[ii];
    for (ii=0; ii < 7; ii++ ) dprim2[ii] = d7prim2[ii];
    return Rcpp::wrap( 0.1*std::sqrt(6./7.*D7Dist(dprim1,dprim2)));
}


