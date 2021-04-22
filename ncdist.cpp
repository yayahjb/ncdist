
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Mat66.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "NCDist.h"
#include <stdlib.h>
#include "S6M_SellingReduce.h"

int dodc7=0, dodc7sq=0;

//*****************************************************************************
G6 makeprimredprobe( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma )
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
            mc = rawcell.LatSymMat66(latsym);
            std::cout << mc << std::endl;
            primcell = mc*(rawcell.Cell2V6());
            std::cout << primcell << std::endl;
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
    CS6M_G6Reduce(dprimcell,dg6redprimcell,reduced);
    /*    std::cout << "dprimcell: " << dprimcell[0]<<" "<< dprimcell[1]<<" "<< dprimcell[2]<<" "<< dprimcell[3]<<" "<< dprimcell[4]<<" "<< dprimcell[5]<<" " << std::endl;
    std::cout << "dg6redprimcell: " << dg6redprimcell[0]<<" "<< dg6redprimcell[1]<<" "<< dg6redprimcell[2]<<" "
       << dg6redprimcell[3]<<" "<< dg6redprimcell[4]<<" "<< dg6redprimcell[5]<<" " << std::endl;
    */
    g6redprimcell = G6(dg6redprimcell);
    return g6redprimcell;
}


int main(int argc, char ** argv) {

    std::string lat1, lat2;
    double a1,b1,c1,alpha1,beta1,gamma1;
    double a2,b2,c2,alpha2,beta2,gamma2;
    G6 prim1, prim2;
    double dprim1[6];
    double dprim2[6];
    double dc71[7],dc73[7];
    double dc72[7],dc74[7];
    double dcdist;
    size_t ii;
    dodc7 = dodc7sq = 0;
    if (argc < 15) {
        std::cerr 
		<< "Usage: ncdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--dodc7]" 
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
    if (argc > 15 && std::string(argv[15]) == "--dodc7") dodc7=1;
    if (argc > 15 && std::string(argv[15]) == "--dodc7sq") dodc7sq=1;
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    /* std::cout << "prim1: "<<prim1 << std::endl;
    std::cout << "prim2: "<<prim2 << std::endl; */
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    CS6M_G6toDC7(dprim1,dc71);
    CS6M_G6toDC7(dprim2,dc72);
    dcdist=0.;
    for (ii=0; ii<7; ii++) {
      if (dodc7) {
        dc73[ii]=std::sqrt(dc71[ii]);
        dc74[ii]=std::sqrt(dc72[ii]);
        dcdist+=(dc73[ii]-dc74[ii])*(dc73[ii]-dc74[ii]);
      } else if (dodc7sq) {
        dcdist+=(dc71[ii]-dc72[ii])*(dc71[ii]-dc72[ii]);
      }
    }
    std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "
              << dprim1[3] << ", " << dprim1[4] << ", " << dprim1[5] <<"]" << std::endl;
    std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "
              << dprim2[3] << ", " << dprim2[4] << ", " << dprim2[5] <<"]" << std::endl;
    std::cout << 0.1*std::sqrt(NCDist(dprim1,dprim2)) << std::endl;
    if (dodc7) {
      std::cout << "dc73: [" << dc73[0] <<", "<< dc73[1] << ", "<< dc73[2] << ", "
              << dc73[3] << ", " << dc73[4] << ", " << dc73[5] << ", " << dc73[6] << "]" << std::endl;
      std::cout << "dc74: [" << dc74[0] <<", "<< dc74[1] << ", "<< dc74[2] << ", "
              << dc74[3] << ", " << dc74[4] << ", " << dc74[5] << ", " << dc74[6] <<"]" << std::endl;
      std::cout << std::sqrt(dcdist/7) << std::endl;
    }
    if (dodc7sq) {
      std::cout << "dc71: [" << dc71[0] <<", "<< dc71[1] << ", "<< dc71[2] << ", "
              << dc71[3] << ", " << dc71[4] << ", " << dc71[5] << ", " << dc71[6] << "]" << std::endl;
      std::cout << "dc72: [" << dc72[0] <<", "<< dc72[1] << ", "<< dc72[2] << ", "
              << dc72[3] << ", " << dc72[4] << ", " << dc72[5] << ", " << dc72[6] <<"]" << std::endl;
      std::cout << 0.1*std::sqrt(std::sqrt(dcdist)) << std::endl;
    }
    return 0;
}
