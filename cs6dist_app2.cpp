
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Reducer.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
#include "CS6Dist_func.h"
#include "S6M_SellingReduce.h"
#include "S6_primredcell.h"

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
    G6 dredprimcell;
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
           primcell[5] = extra-a-b;
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
            std::cerr << "Unrecognized lattice symbol "<< testlattice<<" treated as P" << std::endl;
            latsym = "P";
            CS6M_CelltoG6(rawcell,v6cell);
            CS6M_LatSymMat66(v6cell,clatsym,mc,primcell);
            /*mc = rawcell.LatSymMat66(latsym);
            primcell = mc*(rawcell.Cell2V6());*/
            break;
    }
    reduced=0;
    CS6M_G6Reduce(primcell,redprimcell,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) redprimcell[ii]=redprimcell_as_g6[ii]=primcell[ii];
    } else {
      for(ii=0;ii<6;ii++) redprimcell_as_g6[ii]=redprimcell[ii];
    }
    CS6M_G6toD7(primcell,d7primcell);
    reduced=0;
    CS6M_D7Reduce(d7primcell,d7redprimcell,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) d7redprimcell_as_g6[ii]=primcell[ii];
      for(ii=0;ii<7;ii++) d7redprimcell[ii]=d7primcell[ii];
    } else {
      CS6M_D7toG6(d7redprimcell,d7redprimcell_as_g6);
    }
    CS6M_G6toS6(primcell,s6primcell);
    reduced=0;
    CS6M_S6Reduce(s6primcell,s6redprimcell,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) s6redprimcell[ii]=s6redprimcell_as_g6[ii]=s6primcell[ii];
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
    double a2,b2,c2,alpha2,beta2,gamma2,extra2;
    double cellparams1[6], cellparams2[6];
    double * pcellparams1=cellparams1;
    double * pcellparams2=cellparams2;
    S6 prim1, prim2;
    double dprim1[6];
    double * pdprim1;
    double * qdprim1;
    double dprim2[6];
    double * pdprim2;
    double * qdprim2;
    double g6primcell1[6],g6primcell2[6];
    double * pg6primcell1 = g6primcell1;
    double * pg6primcell2 = g6primcell2;
    double s6primcell1[6],s6primcell2[6];
    double * ps6primcell1 = s6primcell1;
    double * ps6primcell2 = s6primcell2;
    double Mat661[36],Mat662[36];
    double * pMat661 = Mat661;
    double * pMat662 = Mat662;
    double celldeg[6];
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

    cellparams1[0] = a1 = atof(argv[2]);
    cellparams2[0] = a2 = atof(argv[9+argoff]);
    cellparams1[1] = b1 = atof(argv[3]);
    cellparams2[1] = b2 = atof(argv[10+argoff]);
    cellparams1[2] = c1 = atof(argv[4]);
    cellparams2[2] = c2 = atof(argv[11+argoff]);
    cellparams1[3] = alpha1 = atof(argv[5]);
    cellparams2[3] = alpha2 = atof(argv[12+argoff]);
    cellparams1[4] = beta1 = atof(argv[6]);
    cellparams2[4] = beta2 = atof(argv[13+argoff]);
    cellparams1[5] = gamma1 = atof(argv[7]);
    cellparams2[5] = gamma2 = atof(argv[14+argoff]);
    
    prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1);
    prim2 = makeprimredcell(lat2,a2,b2,c2,alpha2,beta2,gamma2,extra2);
    pdprim1 = dprim1;
    pdprim2 = dprim2;
    qdprim1 = s6_primredcell(&(lat1[0]),pcellparams1,
        0,pg6primcell1,ps6primcell1,pdprim1,Mat661);
    qdprim2 = s6_primredcell(&(lat2[0]),pcellparams2,
        0,pg6primcell2,ps6primcell2,pdprim2,Mat662);
    

    if (!qdprim1) std::cout << "qdprim1 is false" << std::endl;
    if (!qdprim2) std::cout << "qdprim2 is false" << std::endl;

    for (ii=0; ii < 6; ii++) {
      if (fabs(dprim1[ii] - prim1[ii])> 1.e-10) 
        std::cout << "prim1 mismatch " << ii << " " 
        << dprim1[ii] << "  " << prim1[ii] << " " << dprim1[ii]-prim1[ii] << std::endl;
      if (fabs(dprim2[ii] - prim2[ii]) > 1.e-10) 
        std::cout << "prim2 mismatch " << ii << " " 
        << dprim2[ii] << "  " << prim2[ii] << " " << dprim2[ii]-prim2[ii] << std::endl;
    }


    std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "
              << dprim1[3] << ", " << dprim1[4] << ", " << dprim1[5] <<  "]" << std::endl;
    std::cout << "g6primcell1: [" << g6primcell1[0] <<", "<< g6primcell1[1] << ", "<< g6primcell1[2] << ", "
              << g6primcell1[3] << ", " << g6primcell1[4] << ", " << g6primcell1[5] <<  "]" << std::endl;
    std::cout << "s6primcell1: [" << s6primcell1[0] <<", "<< s6primcell1[1] << ", "<< s6primcell1[2] << ", "
              << s6primcell1[3] << ", " << s6primcell1[4] << ", " << s6primcell1[5] <<  "]" << std::endl;
    CS6M_G6toCelldeg(g6primcell1,celldeg);
    std::cout << "g6primcelldeg1: ["<< celldeg[0] << ", "  << celldeg[1] << ", "
        << celldeg[2] << ", "<< celldeg[3] << ", "<< celldeg[4] << ", "<< celldeg[5] << "]" << std::endl;
    CS6M_S6toG6(s6primcell1,g6primcell1);
    CS6M_G6toCelldeg(g6primcell1,celldeg);
    std::cout << "s6primcelldeg1: ["<< celldeg[0] << ", "  << celldeg[1] << ", "
        << celldeg[2] << ", "<< celldeg[3] << ", "<< celldeg[4] << ", "<< celldeg[5] << "]" << std::endl;
    std::cout << "Mat661 ["<< Mat661[0]  << ", "<< Mat661[1]  << ", "<< Mat661[2]  << ", "<< Mat661[3]  << ", "<< Mat661[4]  << ", "<< Mat661[5]  << "," << std::endl;
    std::cout << "        "<< Mat661[6]  << ", "<< Mat661[7]  << ", "<< Mat661[8]  << ", "<< Mat661[9]  << ", "<< Mat661[10] << ", "<< Mat661[11] << "," << std::endl;
    std::cout << "        "<< Mat661[12] << ", "<< Mat661[13] << ", "<< Mat661[14] << ", "<< Mat661[15] << ", "<< Mat661[16] << ", "<< Mat661[17] << "," << std::endl;
    std::cout << "        "<< Mat661[18] << ", "<< Mat661[19] << ", "<< Mat661[20] << ", "<< Mat661[21] << ", "<< Mat661[22] << ", "<< Mat661[23] << "," << std::endl;
    std::cout << "        "<< Mat661[24] << ", "<< Mat661[25] << ", "<< Mat661[26] << ", "<< Mat661[27] << ", "<< Mat661[28] << ", "<< Mat661[28] << "," << std::endl;
    std::cout << "        "<< Mat661[30] << ", "<< Mat661[31] << ", "<< Mat661[32] << ", "<< Mat661[33] << ", "<< Mat661[34] << ", "<< Mat661[35] << "]" << std::endl;
    std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "
              << dprim2[3] << ", " << dprim2[4] << ", " << dprim2[5] <<  "]" << std::endl;
    std::cout << "g6primcell2: [" << g6primcell2[0] <<", "<< g6primcell2[1] << ", "<< g6primcell2[2] << ", "
              << g6primcell2[3] << ", " << g6primcell2[4] << ", " << g6primcell2[5] <<  "]" << std::endl;
    std::cout << "s6primcell2: [" << s6primcell2[0] <<", "<< s6primcell2[1] << ", "<< s6primcell1[2] << ", "
              << s6primcell2[3] << ", " << s6primcell2[4] << ", " << s6primcell2[5] <<  "]" << std::endl;
    CS6M_G6toCelldeg(g6primcell2,celldeg);
    std::cout << "g6primcelldeg2: ["<< celldeg[0] << ", "  << celldeg[1] << ", "
        << celldeg[2] << ", "<< celldeg[3] << ", "<< celldeg[4] << ", "<< celldeg[5] << "]" << std::endl;
    CS6M_S6toG6(s6primcell2,g6primcell2);
    CS6M_G6toCelldeg(g6primcell2,celldeg);
    std::cout << "s6primcelldeg2: ["<< celldeg[0] << ", "  << celldeg[1] << ", "
        << celldeg[2] << ", "<< celldeg[3] << ", "<< celldeg[4] << ", "<< celldeg[5] << "]" << std::endl;
    std::cout << "Mat662 ["<< Mat662[0]  << ", "<< Mat662[1]  << ", "<< Mat662[2]  << ", "<< Mat662[3]  << ", "<< Mat662[4]  << ", "<< Mat662[5]  << "," << std::endl;
    std::cout << "        "<< Mat662[6]  << ", "<< Mat662[7]  << ", "<< Mat662[8]  << ", "<< Mat662[9]  << ", "<< Mat662[10] << ", "<< Mat662[11] << "," << std::endl;
    std::cout << "        "<< Mat662[12] << ", "<< Mat662[13] << ", "<< Mat662[14] << ", "<< Mat662[15] << ", "<< Mat662[16] << ", "<< Mat662[17] << "," << std::endl;
    std::cout << "        "<< Mat662[18] << ", "<< Mat662[19] << ", "<< Mat662[20] << ", "<< Mat662[21] << ", "<< Mat662[22] << ", "<< Mat662[23] << "," << std::endl;
    std::cout << "        "<< Mat662[24] << ", "<< Mat662[25] << ", "<< Mat662[26] << ", "<< Mat662[27] << ", "<< Mat662[28] << ", "<< Mat662[28] << "," << std::endl;
    std::cout << "        "<< Mat662[30] << ", "<< Mat662[31] << ", "<< Mat662[32] << ", "<< Mat662[33] << ", "<< Mat662[34] << ", "<< Mat662[35] << "]" << std::endl;
    rawdist = CS6Dist_func(dprim1,dprim2);
    std::cout << "rawdist: "<< rawdist << std::endl;
    std::cout << 0.1*std::sqrt(rawdist) << std::endl;
    return 0;
}
