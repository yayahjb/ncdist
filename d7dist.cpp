

#include "S6M_SellingReduce.h"
#include "D7Dist.h"
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
#include "Delone.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <ctype.h>

//*****************************************************************************
G6 makeprimredprobe( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma )
{
    int reduced;
    std::string latsym;
    char clatsym;
    G6 v6cell;
    G6 redprimcell;
    G6 dredprimcell;
    Mat66 mc;
    Mat66 m;
    Mat66 dm;
    G6 primcell;
    G6 recipcell;
    G6 reducedBase;
    G6 primredprobe;
    G6 dprimredprobe;
    double d7cell[7];
    double d7redcell[7];
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
        default:
            std::cerr << "Unrecognized lattice symbol "<< testlattice<<" treated as P" << std::endl;
            latsym = "P";
            clatsym = 'P';
            CS6M_CelltoG6(rawcell,v6cell);
            CS6M_LatSymMat66(v6cell,clatsym,mc,primcell);
            break;
    }
    CS6M_G6Reduce(primcell,redprimcell,reduced);
    CS6M_G6toD7(primcell,d7cell);
    CS6M_D7Reduce(d7cell,d7redcell,reduced);
    CS6M_D7toG6(d7redcell,dredprimcell);
    return dredprimcell;
}


int main(int argc, char ** argv) {

    std::string lat1, lat2;
    double a1,b1,c1,alpha1,beta1,gamma1;
    double a2,b2,c2,alpha2,beta2,gamma2;
    G6 prim1, prim2;
    double dprim1[7];
    double dprim2[7];
    size_t ii;
         
    if (argc < 15) {
        std::cerr 
		<< "Usage: ncdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2" 
		<< std::endl;
        return -1;
    }
    lat1 = std::string(argv[1]);
    lat2 = std::string(argv[8]);
    a1 = atof(argv[2]);
    a2 = atof(argv[9]);
    b1 = atof(argv[3]);
    b2 = atof(argv[10]);
    c1 = atof(argv[4]);
    c2 = atof(argv[11]);
    alpha1 = atof(argv[5]);
    alpha2 = atof(argv[12]);
    beta1 = atof(argv[6]);
    beta2 = atof(argv[13]);
    gamma1 = atof(argv[7]);
    gamma2 = atof(argv[14]);
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    for (ii=0; ii < 3; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    dprim1[3] = prim1[0]+prim1[1]+prim1[2]+prim1[3]+prim1[4]+prim1[5];
    dprim2[3] = prim2[0]+prim2[1]+prim2[2]+prim2[3]+prim2[4]+prim2[5];
    dprim1[4] =         +prim1[1]+prim1[2]+prim1[3];
    dprim2[4] =         +prim2[1]+prim2[2]+prim2[3];
    dprim1[5] = prim1[0]         +prim1[2]         +prim1[4];
    dprim2[5] = prim2[0]         +prim2[2]         +prim2[4];
    dprim1[6] = prim1[0]+prim1[1]                           +prim1[5];
    dprim2[6] = prim2[0]+prim2[1]                           +prim2[5];
    std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "
              << dprim1[3] << ", " << dprim1[4] << ", " << dprim1[5] << ", " << dprim1[6]<< "]" << std::endl;
    std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "
              << dprim2[3] << ", " << dprim2[4] << ", " << dprim2[5] << ", " << dprim1[6]<< "]" << std::endl;
    std::cout << 0.1*std::sqrt(D7Dist(dprim1,dprim2)) << std::endl;
    return 0;
}
