#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Mat66.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "NCDist.h"
#include <stdlib.h>
#define CS6M_DEBUG zz
#include "S6M_SellingReduce.h"

int donc=0, dodc7=0, dodc7sq=0, dodc7unsrt=0, dodc7unsrtsq=0, dodc10=0, dodc10sq=0, info=0;

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
    double dc13[13];
    double dc7[7];
    double dc7_unsrt[7];
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
    if (info) {
      std::cout << "dprimcell: " << dprimcell[0]<<" "<< dprimcell[1]<<" "<< dprimcell[2]<<" "<< dprimcell[3]<<" "<< dprimcell[4]<<" "<< dprimcell[5]<<" " << std::endl;
      std::cout << "dg6redprimcell: " << dg6redprimcell[0]<<" "<< dg6redprimcell[1]<<" "<< dg6redprimcell[2]<<" "
       << dg6redprimcell[3]<<" "<< dg6redprimcell[4]<<" "<< dg6redprimcell[5]<<" " << std::endl;
      CS6M_G6toDC13blocks(dg6redprimcell,dc13);
      std::cout << "dc13_blocks: " << std::sqrt(dc13[0]) << " "<< std::sqrt(dc13[1]) << " "<< std::sqrt(dc13[2]) << " "<< std::sqrt(dc13[3]) 
                                   << " "<< std::sqrt(dc13[4]) << " "<< std::sqrt(dc13[5]) << " "<< std::sqrt(dc13[6]) << " "<< std::sqrt(dc13[7])
                                   << " "<< std::sqrt(dc13[8]) << " "<< std::sqrt(dc13[9]) << " "<< std::sqrt(dc13[10]) << " "<< std::sqrt(dc13[11])
                                   << " "<< std::sqrt(dc13[12]) << " " << std::endl;
      CS6M_G6toDC7(dg6redprimcell,dc7);
      std::cout << "dc7: " << std::sqrt(dc7[0]) << " "<< std::sqrt(dc7[1]) << " "<< std::sqrt(dc7[2]) << " "<< std::sqrt(dc7[3]) 
                     << " "<< std::sqrt(dc7[4]) << " "<< std::sqrt(dc7[5]) << " "<< std::sqrt(dc7[6]) << " " << std::endl;
      CS6M_G6toDC7unsorted(dg6redprimcell,dc7_unsrt);
      std::cout << "dc7_unsorted: " <<  std::sqrt(dc7_unsrt[0]) << " " <<  std::sqrt(dc7_unsrt[1]) << " " <<  std::sqrt(dc7_unsrt[2]) << " " <<  std::sqrt(dc7_unsrt[3])
                                    << " " <<  std::sqrt(dc7_unsrt[4]) << " " <<  std::sqrt(dc7_unsrt[5]) << " " <<  std::sqrt(dc7_unsrt[6]) << " " << std::endl;
    }
    g6redprimcell = G6(dg6redprimcell);
    return g6redprimcell;
}


int main(int argc, char ** argv) {

    std::string lat1, lat2;
    double a1,b1,c1,alpha1,beta1,gamma1;
    double a2,b2,c2,alpha2,beta2,gamma2;
    G6 prim1, prim2;
    double delta1, delta2;
    double dprim1[6];
    double dprim2[6];
    size_t ii;
    std::size_t notfound=std::string::npos;
    dodc7 = dodc7sq = dodc7unsrt = dodc7unsrtsq = dodc10 = dodc10sq = donc = info = 0;
    if (std::string(argv[0]).rfind("ncdist") != notfound) {
      donc = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: ncdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 \\" << std::endl
                << "[--info] [--dodc7] [--dodc7sq] [--dodc7unsrt] [--dodc7unsrtsq] [--dodc10] [--dodc10sq]" 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc7dist") != notfound) {
      dodc7 = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc7dist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc7sq]" 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc7sqdist") != notfound ) {
      dodc7sq = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc7sqdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc7]" 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc7unsrtdist") != notfound ) {
      dodc7sq = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc7unsrtdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc7unsrtsq]" 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc7unsrtsqdist") != notfound ) {
      dodc7sq = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc7unsrtsqdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc7unsrt]" 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc10dist") != notfound) {
      dodc10 = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc10dist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc10sq] " 
  		<< std::endl;
          return -1;
      }
    }
    if (std::string(argv[0]).rfind("dc10sqdist") != notfound ) {
      dodc10sq = 1;
      if (argc < 15) {
          std::cerr 
  		<< "Usage: dc10sqdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2 [--info] [--donc] [--dodc10]" 
  		<< std::endl;
          return -1;
      }
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
    for (ii=15; ii < argc; ii++) {
      if (std::string(argv[ii]) == "--info") info=1;
      if (std::string(argv[ii]) == "--donc") donc=1;
      if (std::string(argv[ii]) == "--dodc7") dodc7=1;
      if (std::string(argv[ii]) == "--dodc7sq") dodc7sq=1;
      if (std::string(argv[ii]) == "--dodc7unsrt") dodc7unsrt=1;
      if (std::string(argv[ii]) == "--dodc7unsrtsq") dodc7unsrtsq=1;
      if (std::string(argv[ii]) == "--dodc10") dodc10=1;
      if (std::string(argv[ii]) == "--dodc10sq") dodc10sq=1;
    }
    prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1);
    prim2 = makeprimredprobe(lat2,a2,b2,c2,alpha2,beta2,gamma2);
    for (ii=0; ii < 6; ii++) {
      dprim1[ii] = prim1[ii];
      dprim2[ii] = prim2[ii];
    }
    delta1 = dprim1[0]*1.e-9;
    delta2 = dprim2[0]*1.e-9;
    if (info) {
      std::cout << "dprim1: [" << dprim1[0] <<", "<< dprim1[1] << ", "<< dprim1[2] << ", "
              << dprim1[3] << ", " << dprim1[4] << ", " << dprim1[5] <<"]" << std::endl;
      std::cout << "dprim2: [" << dprim2[0] <<", "<< dprim2[1] << ", "<< dprim2[2] << ", "
              << dprim2[3] << ", " << dprim2[4] << ", " << dprim2[5] <<"]" << std::endl;
    }
    if (donc && info) std::cout << "raw ncdist: " << NCDist(dprim1,dprim2) << std::endl;
    if (donc) std::cout << 0.1*std::sqrt(NCDist(dprim1,dprim2)) << std::endl;
    if (dodc7 && info) std::cout << "raw dc7 dist: "<< DC7Distraw(dprim1,dprim2) << std::endl;
    if (dodc7) std::cout << DC7Dist(dprim1,dprim2) << std::endl;
    if (dodc7sq && info) std::cout << "raw dc7sq dist: "<<DC7sqDistraw(dprim1,dprim2) << std::endl;
    if (dodc7sq) std::cout << DC7sqDist(dprim1,dprim2) << std::endl;
    if (dodc7unsrt && info) std::cout << "raw dc7unsrt dist: "<< DC7unsrtDistraw(dprim1,dprim2) << std::endl;
    if (dodc7unsrt) std::cout << DC7unsrtDist(dprim1,dprim2) << std::endl;
    if (dodc7unsrtsq && info) {
      double dc71[7], dc72[7];
      double g6vec1[6], g6vec2[6];
      int error1, error2;
      CS6M_G6toDC7unsorted(dprim1,dc71);
      CS6M_G6toDC7unsorted(dprim2,dc72);
      CS6M_DC7unsortedtoG6(dc71,g6vec1,error1);
      CS6M_DC7unsortedtoG6(dc72,g6vec2,error2);
      if (!error1 &&
          (CS6M_abs(g6vec1[0]-dprim1[0])>delta1 
           || CS6M_abs(g6vec1[1]-dprim1[1])>delta1
           || CS6M_abs(g6vec1[2]-dprim1[2])>delta1
           || CS6M_abs(g6vec1[3]-dprim1[3])>delta1
           || CS6M_abs(g6vec1[4]-dprim1[4])>delta1
           || CS6M_abs(g6vec1[5]-dprim1[5])>delta1)
         ) {
        std::cout << "g6vec1 ["<< g6vec1[0] << " " << g6vec1[1] << " " << g6vec1[2] << " " << g6vec1[3] << " " << g6vec1[4] << " " << g6vec1[5] << " "<< std::endl;
      }
      if (!error2 &&   
          (CS6M_abs(g6vec2[0]-dprim2[0])>delta1      
           || CS6M_abs(g6vec2[1]-dprim2[1])>delta1
           || CS6M_abs(g6vec2[2]-dprim2[2])>delta1
           || CS6M_abs(g6vec2[3]-dprim2[3])>delta1
           || CS6M_abs(g6vec2[4]-dprim2[4])>delta1
           || CS6M_abs(g6vec2[5]-dprim2[5])>delta1)
         ) {
        std::cout << "g6vec2 ["<< g6vec2[0] << " " << g6vec2[1] << " " << g6vec2[2] << " " << g6vec2[3] << " " << g6vec2[4] << " " << g6vec2[5] << " "<< std::endl;
      }
      std::cout << "raw dc7unsrtsq dist: "<<DC7unsrtsqDistraw(dprim1,dprim2) << std::endl;
    }
    if (dodc7unsrtsq) std::cout << DC7unsrtsqDist(dprim1,dprim2) << std::endl;
    if (dodc10 && info) std::cout << "raw dc10 dist: "<< DC10Distraw(dprim1,dprim2) << std::endl;
    if (dodc10) std::cout << DC10Dist(dprim1,dprim2) << std::endl;
    if (dodc10sq && info) std::cout << "raw dc10sq dist: "<<DC10sqDistraw(dprim1,dprim2) << std::endl;
    if (dodc10sq && info) {
      double dc131[13], dc132[13];
      double g6vecs1[6][6], g6vecs2[6][6];
      int error1,error2,ng6vecs1,ng6vecs2;
      int ii;
      CS6M_G6toDC13blocks(dprim1,dc131);
      CS6M_G6toDC13blocks(dprim2,dc132);
      CS6M_DC13blockstoG6(dc131,g6vecs1,6,ng6vecs1,error1);
      CS6M_DC13blockstoG6(dc132,g6vecs2,6,ng6vecs2,error2);
      if (!error1 && ng6vecs1 > 0) {
        for (ii=0; ii < ng6vecs1; ii++) {
          std::cout << "g6vecs1["<<ii<<"] ["<<g6vecs1[0][ii]<<","<<g6vecs1[1][ii]<<","<<g6vecs1[2][ii]<<","
            <<g6vecs1[3][ii]<<","<<g6vecs1[4][ii]<<","<<g6vecs1[5][ii]<<"]"<<std::endl;
        }
      } else {
        std::cout <<  "CS6M_DC13blockstoG6(dc131,g6vecs1,720,ng6vecs1,error1) failed, error1= "<< error1 << ", ng6vecs1= "<< ng6vecs1<<std::endl;
      }
      if (!error2 && ng6vecs2 > 0) {
        for (ii=0; ii < ng6vecs2; ii++) {
          std::cout << "g6vecs2["<<ii<<"] ["<<g6vecs2[0][ii]<<","<<g6vecs2[1][ii]<<","<<g6vecs2[2][ii]<<","
            <<g6vecs2[3][ii]<<","<<g6vecs2[4][ii]<<","<<g6vecs2[5][ii]<<"]"<<std::endl;
        }
      } else {
        std::cout <<  "CS6M_DC13blockstoG6(dc132,g6vecs2,720,ng6vecs2,error2) failed, error2= "<< error2 << ", ng6vecs2= "<< ng6vecs2<<std::endl;
      }
    }
    if (dodc10sq) std::cout << DC10sqDist(dprim1,dprim2) << std::endl;
    return 0;
}
