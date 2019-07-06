
#include <RcppParallel.h>
using namespace RcppParallel;
#include <RcppArmadillo.h>

#include "S6M_SellingReduce.h"
#include "Reducer.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h" 
#include "NCDist.h"
#include <stdlib.h>

//*****************************************************************************
G6 makeprimredprobe( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma )
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 redprimcell;
    double dredprimcell[6];
    Mat66 mc;
    Mat66 m;
    G6 primcell;
    double dprimcell[6];
    double dreduced;
    G6 recipcell;
    G6 reducedBase;
    G6 primredprobe;
    double crootvol;
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
            mc = rawcell.LatSymMat66(latsym);
            primcell = mc*(rawcell.Cell2V6());
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
        default:
            /* Rprintf("Unrecognized lattice symbol %s treated as P\n",testlattice.c_str()); */
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
    CS6M_G6Reduce(dprimcell,dredprimcell,dreduced);
    redprimcell = G6(dredprimcell);
    return redprimcell;
}

extern "C" SEXP rcpp_ncdist ( SEXP lat1_, SEXP a1_, SEXP b1_, SEXP c1_, 
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
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    //std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "<< dprim1[3] << ", " << dprim1[4] << ", " << dprim1[1] <<"]" << std::endl;
    //std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "<< dprim2[3] << ", " << dprim2[4] << ", " << dprim2[1] <<"]" << std::endl;
    return Rcpp::wrap( 0.1*std::sqrt(NCDist(dprim1,dprim2)));
}


