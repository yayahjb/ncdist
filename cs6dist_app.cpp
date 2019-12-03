
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Reducer.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
#include "CS6Dist_func.h"
#include "S6M_SellingReduce.h"

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <ctype.h>

// *****************************************************************************
S6 makeprimredcell( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma, double extra )
{
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 redprimcell;
    G6 redprimcell_as_g6;
    D7 d7redprimcell;
    G6 d7redprimcell_as_g6;
    S6 s6redprimcell;
    G6 s6redprimcell_as_g6;
    double d7primcell[7];
    double s6primcell[6];
    double dredprimcell[6];
    Mat66 mc;
    double dmc[36];
    Mat66 m;
    Mat66 dm;
    G6 primcell;
    double dprimcell[6];
    G6 recipcell;
    G6 reducedBase;
    G6 g6primredprobe_as_g6;
    D7 d7primredprobe_as_g6;
    S6 S6primredprobe_as_g6;
    LRL_Cell g6primredprobe; 
    LRL_Cell d7primredprobe; 
    LRL_Cell s6primredprobe; 
    double crootvol;
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
            CS6M_LatSymMat66(v6cell,clatsym,dmc,dprimcell);
            mc=Mat66(dmc);
            primcell=G6(dprimcell);
            break;
        case 'V':
        case 'v':
	    dprimcell[0] = a;
	    dprimcell[1] = b;
	    dprimcell[2] = c;
	    dprimcell[3] = alpha;
	    dprimcell[4] = beta;
	    dprimcell[5] = gamma;
            primcell=G6(dprimcell);
            break;
        case 'D':
        case 'd':
           dprimcell[0] = a;
           dprimcell[1] = b;
           dprimcell[2] = c;
           dprimcell[3] = beta-b-c;
           dprimcell[4] = gamma-a-c;
           dprimcell[5] = extra-a-b;
           primcell=G6(dprimcell);
           break;
        case 'S':
        case 's':
           dprimcell[3] = 2.*a;
           dprimcell[4] = 2.*b;
           dprimcell[5] = 2.*c;
           dprimcell[0] = -alpha-c-b;
           dprimcell[1] = -beta-c-a;
           dprimcell[2] = -gamma-b-a;
           primcell=G6(dprimcell);
           break;
        default:
            std::cerr << "Unrecognized lattice symbol "<< testlattice<<" treated as P" << std::endl;
            latsym = "P";
            CS6M_CelltoG6(rawcell,v6cell);
            CS6M_LatSymMat66(v6cell,clatsym,dmc,dprimcell);
            mc=Mat66(dmc);
            primcell=G6(dprimcell);
            break;
    }
    reduced=0;
    CS6M_G6Reduce(dprimcell,dredprimcell,reduced);
    redprimcell=G6(dredprimcell);
    if (!reduced) {
      for(ii=0;ii<6;ii++) redprimcell[ii]=redprimcell_as_g6[ii]=0.;
    } else {
      for(ii=0;ii<6;ii++) redprimcell_as_g6[ii]=redprimcell[ii];
    }
    CS6M_G6toD7(primcell,d7primcell);
    reduced=0;
    CS6M_D7Reduce(d7primcell,d7redprimcell,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) d7redprimcell_as_g6[ii]=0.;
      for(ii=0;ii<7;ii++) d7redprimcell[ii]=0.;
    } else {
      CS6M_D7toG6(d7redprimcell,d7redprimcell_as_g6);
    }
    CS6M_G6toS6(primcell,s6primcell);
    reduced=0;
    CS6M_S6Reduce(s6primcell,s6redprimcell,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) s6redprimcell[ii]=s6redprimcell_as_g6[ii]=0;
    } else {
      CS6M_S6toG6(s6redprimcell,s6redprimcell_as_g6);
    }
    g6primredprobe = LRL_Cell_Degrees(redprimcell_as_g6);
    d7primredprobe = LRL_Cell_Degrees(d7redprimcell_as_g6);
    s6primredprobe = LRL_Cell_Degrees(s6redprimcell_as_g6);
    std::cout << "Primitive Reduced Probe LRL_Cell: " <<
    g6primredprobe[0]<<" "<<
    g6primredprobe[1]<<" "<<
    g6primredprobe[2]<<" "<<
    g6primredprobe[3]<<" "<<
    g6primredprobe[4]<<" "<<
    g6primredprobe[5] << std::endl;
    std::cout << "Delaunay Primitive Reduced Probe LRL_Cell: " <<
    d7primredprobe[0]<<" "<<
    d7primredprobe[1]<<" "<<
    d7primredprobe[2]<<" "<<
    d7primredprobe[3]<<" "<<
    d7primredprobe[4]<<" "<<
    d7primredprobe[5] << std::endl;
    std::cout << "Selling Primitive Reduced Probe S6 LRL_Cell: " <<
    s6redprimcell[0]<<" "<<
    s6redprimcell[1]<<" "<<
    s6redprimcell[2]<<" "<<
    s6redprimcell[3]<<" "<<
    s6redprimcell[4]<<" "<<
    s6redprimcell[5] << std::endl;
    std::cout << "Volume :" << LRL_Cell(redprimcell).Volume() << std::endl;
    crootvol = pow(LRL_Cell(redprimcell).Volume(),1./3.);
    Reducer::Reduce((LRL_Cell(redprimcell).Inverse()).Cell2V6(),m,reducedBase,0.0);
    recipcell = LRL_Cell_Degrees((LRL_Cell(redprimcell).Inverse()));
    std::cout << "Reciprocal of Primitive Probe LRL_Cell: " <<
    recipcell[0]<<" "<<
    recipcell[1]<<" "<<
    recipcell[2]<<" "<<
    recipcell[3]<<" "<<
    recipcell[4]<<" "<<
    recipcell[5]<< std::endl;
    std::cout << "Volume of Reciprocal LRL_Cell: " <<
    (LRL_Cell(redprimcell).Inverse()).Volume() << std::endl;
    std::cout << "V7 linearized and scaled: "
    << g6primredprobe[0]*std::sqrt(6./7.)<<" "
    << g6primredprobe[1]*std::sqrt(6./7.)<<" "
    << g6primredprobe[2]*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[0])*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[1])*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[2])*std::sqrt(6./7.)<<" "
    << " "<<
    crootvol*std::sqrt(6./7.)<<std::endl;
    if (latsym[0] == 'V' || latsym[0] == 'v') {
        std::cout << "raw G6 vector: "
        << primcell[0]<<" "
        << primcell[1]<<" "
        << primcell[2]<<" "
        << primcell[3]<<" "
        << primcell[4]<<" "
        << primcell[5]<<std::endl;
    } else {
        std::cout << "raw G6 vector for S6: "
        << s6primredprobe[0]*s6primredprobe[0]<<" "
        << s6primredprobe[1]*s6primredprobe[1]<<" "
        << s6primredprobe[2]*s6primredprobe[2]<<" "
        << 2.*s6primredprobe[1]*s6primredprobe[2]*cos(s6primredprobe[3]*std::atan(1.0)/45.)<<" "
        << 2.*s6primredprobe[0]*s6primredprobe[2]*cos(s6primredprobe[4]*std::atan(1.0)/45.)<<" "
        << 2.*s6primredprobe[0]*s6primredprobe[1]*cos(s6primredprobe[5]*std::atan(1.0)/45.)<<std::endl;
    }
    std::cout << std::endl;
    return S6(s6redprimcell);
}


int main(int argc, char ** argv) {

    std::string lat1, lat2;
    double a1,b1,c1,alpha1,beta1,gamma1,extra1;
    double a2,b2,c2,alpha2,beta2,gamma2,extra2;;
    S6 prim1, prim2;
    double dprim1[6];
    double dprim2[6];
    double rawdist;
    size_t ii;
    char clatsym;
    int argoff;

         
    if (argc < 15) {
        std::cerr 
		<< "Usage: cs6dist_app lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2" 
		<< std::endl;
        return -1;
    }
    clatsym=((std::string(argv[8])).substr(0,1))[0];
    argoff = 1;
    if (isalpha(clatsym)) argoff = 0;
    extra1 = 0.;
    extra2 = 0.;
    lat1 = std::string(argv[1]);
    clatsym= lat1.substr(0,1)[0];
    /* std::cerr <<  "clatsym, argoff " << clatsym << " " << argoff <<  std::endl; */
    if (clatsym == 'D' || clatsym == 'd') {
        argoff = 1;
        if (argc < 16) {
            std::cerr 
                << "Usage: cs6dist_app lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2" 
                << std::endl;
            return -1;
        }
        extra1 = atof(argv[8]);
    }
    lat2 = std::string(argv[8+argoff]);
    clatsym= lat2.substr(0,2)[0];
    if (clatsym == 'D' || clatsym == 'd') {
        if (argc < 16+argoff){
            std::cerr 
                << "Usage: cs6dist_app lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2" 
                << std::endl;
            return -1;
        }
        extra2 = atof(argv[15+argoff]);
    }

    a1 = atof(argv[2]);
    a2 = atof(argv[9+argoff]);
    b1 = atof(argv[3]);
    b2 = atof(argv[10+argoff]);
    c1 = atof(argv[4]);
    c2 = atof(argv[11+argoff]);
    alpha1 = atof(argv[5]);
    alpha2 = atof(argv[12+argoff]);
    beta1 = atof(argv[6]);
    beta2 = atof(argv[13+argoff]);
    gamma1 = atof(argv[7]);
    gamma2 = atof(argv[14+argoff]);
    
    prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1);
    prim2 = makeprimredcell(lat2,a2,b2,c2,alpha2,beta2,gamma2,extra2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "
              << dprim1[3] << ", " << dprim1[4] << ", " << dprim1[5] <<  "]" << std::endl;
    std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "
              << dprim2[3] << ", " << dprim2[4] << ", " << dprim2[5] <<  "]" << std::endl;
    rawdist = CS6Dist_func(dprim1,dprim2);
    std::cout << "rawdist: "<< rawdist << std::endl;
    std::cout << 0.1*std::sqrt(rawdist) << std::endl;
    return 0;
}
