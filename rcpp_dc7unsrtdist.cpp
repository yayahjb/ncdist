
#include <RcppParallel.h>
using namespace RcppParallel;
#include <RcppArmadillo.h> 

#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Mat66.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h" 
#include "NCDist.h"
#include <stdlib.h>
#include "S6M_SellingReduce.h"

//***************************************************************************//
G6 makeprimredprobe( std::string testlattice,
	double a, double b, double c, 
        double alpha, double beta, double gamma, double extra )
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 redprimcell;
    D7 d7redprimcell;
    double dredprimcell[7];
    G6 g6redprimcell;
    double dg6redprimcell[6];
    Mat66 mc;
    Mat66 dm;
    G6 primcell;
    double dprimcell[6];
    double sprimcell[6];
    double dsreduced;
    G6 recipcell;
    G6 reducedBase;
    G6 primredprobe;
    G6 dprimredprobe;
    double crootvol;
    int reduced;
    LRL_Cell rawcell(a,b,c, alpha,beta,gamma);
    int ii;
    bool ret;
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
        case 'V': /* G6 cell */
        case 'v':
	    primcell[0] = a;
	    primcell[1] = b;
	    primcell[2] = c;
	    primcell[3] = alpha;
	    primcell[4] = beta;
	    primcell[5] = gamma;
            break;
        case 'D':  /* D7 cell */
        case 'd':    
           primcell[0] = a;
           primcell[1] = b;
           primcell[2] = c;
           primcell[3] = beta-b-c;
           primcell[4] = gamma-a-c;
           primcell[5] = extra-a-b;
           break;
        case 'S':  /* Selling S6 cell */
        case 's':
           primcell[3] = 2.*a;
           primcell[4] = 2.*b;
           primcell[5] = 2.*c;
           primcell[0] = -alpha-c-b;
           primcell[1] = -beta-c-a;
           primcell[2] = -gamma-b-a;
           break;
        case 'U':  /* unsorted dc7sq cell */
        case 'u':
           primcell[0] = a;
           primcell[1] = b;
           primcell[2] = c;
           {
             double absu, absv, absw, tau;
             absu = b+c-primcell[3];
             absv = a+c-primcell[4];
             absw = a+b-primcell[5];
             if (absu < 0.) absu = 0.;
             if (absv < 0.) absv = 0.;
             if (absw < 0.) absw = 0.;
             tau=a+b+c-absu-absv-absw;
             if (tau < extra) {
               primcell[3] = absu;
               primcell[4] = absv;
               primcell[5] = absw;
             } else {
               primcell[3] = -absu;
               primcell[4] = -absv;
               primcell[5] = -absw;
             }
           }
        default:
            /* Rprintf("Unrecognized lattice symbol %s treated as P\n",testlattice.c_str()) */; 
            latsym = "P";
            mc = rawcell.LatSymMat66(latsym);
            primcell = mc*(rawcell.Cell2V6());
            break;
    }
    dprimcell[0]=primcell[0];
    dprimcell[1]=primcell[1];
    dprimcell[2]=primcell[2];
    dprimcell[3]=primcell[3];
    dprimcell[4]=primcell[4];
    dprimcell[5]=primcell[5];
    reduced=0;
    CS6M_G6Reduce(dprimcell,dg6redprimcell,reduced);
    g6redprimcell = G6(dg6redprimcell);
    /* Rprintf(" dprimcell: [ %g %g %g %g %g %g ]\n", dprimcell[0], dprimcell[1], dprimcell[2], dprimcell[3], dprimcell[4], dprimcell[5]); */
    /* Rprintf("dg6redprimcell: [ %g %g %g %g %g %g ]\n", dg6redprimcell[0], dg6redprimcell[1], dg6redprimcell[2],
       dg6redprimcell[3], dg6redprimcell[4], dg6redprimcell[5]); */

    return g6redprimcell;
}

extern "C" SEXP rcpp_dc7unsrtsqdist ( SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_, 
                              SEXP alpha1_, SEXP beta1_, SEXP gamma1_, SEXP extra1_,
			      SEXP lat2_, SEXP a2_, SEXP b2_, SEXP c2_, 
                              SEXP alpha2_, SEXP beta2_, SEXP gamma2_, SEXP extra2_) {
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
    double extra1    = Rcpp::as<double>(extra1_);
    std::string lat2 = Rcpp::as<std::string>(lat2_);
    double a2        = Rcpp::as<double>(a2_);
    double b2        = Rcpp::as<double>(b2_);
    double c2        = Rcpp::as<double>(c2_);
    double alpha2    = Rcpp::as<double>(alpha2_);
    double beta2     = Rcpp::as<double>(beta2_);
    double gamma2    = Rcpp::as<double>(gamma2_);
    double extra2    = Rcpp::as<double>(extra2_);
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2,extra2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    return Rcpp::wrap(DC7unsrtsqDist(dprim1,dprim2));
}
extern "C" SEXP rcpp_dc7unsrtdist ( SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_, 
                              SEXP alpha1_, SEXP beta1_, SEXP gamma1_, SEXP extra1_,
			      SEXP lat2_, SEXP a2_, SEXP b2_, SEXP c2_, 
                              SEXP alpha2_, SEXP beta2_, SEXP gamma2_, SEXP extra2_) {
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
    double extra1    = Rcpp::as<double>(extra1_);
    std::string lat2 = Rcpp::as<std::string>(lat2_);
    double a2        = Rcpp::as<double>(a2_);
    double b2        = Rcpp::as<double>(b2_);
    double c2        = Rcpp::as<double>(c2_);
    double alpha2    = Rcpp::as<double>(alpha2_);
    double beta2     = Rcpp::as<double>(beta2_);
    double gamma2    = Rcpp::as<double>(gamma2_);
    double extra2    = Rcpp::as<double>(extra2_);
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2,extra2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    return Rcpp::wrap(DC7unsrtDist(dprim1,dprim2));
}


