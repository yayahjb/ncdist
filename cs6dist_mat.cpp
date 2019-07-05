 
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Reducer.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
/* #include "D7Dist.h" */
/* #include "CS6Dist_func.h" */
#include "CS6Dist.h"
#include "S6M_SellingReduce.h"


#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <ctype.h>

using namespace std;

S6 makeprimredcell( std::string testlattice,
        double a, double b, double c, double alpha, double beta, double gamma, double extra, int info=0 )
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
    Mat66 m;
    Mat66 dm;
    G6 primcell;
    G6 recipcell;
    G6 reducedBase;
    G6 g6primredprobe_as_g6;
    D7 d7primredprobe_as_g6;
    S6 S6primredprobe_as_g6;
    LRL_Cell g6primredprobe;
    LRL_Cell d7primredprobe;
    LRL_Cell s6primredprobe;
    g6primredprobe = LRL_Cell_Degrees(redprimcell_as_g6);
    d7primredprobe = LRL_Cell_Degrees(d7redprimcell_as_g6);
    s6primredprobe = LRL_Cell_Degrees(s6redprimcell_as_g6);
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
            CS6M_LatSymMat66(v6cell,clatsym,mc,primcell);
            /* mc = rawcell.LatSymMat66(latsym);
            primcell = mc*(rawcell.Cell2V6());*/
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
    if (info) {
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
   }
    return S6(s6redprimcell);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool not_space( const char c ) {
   return( c != ' ' );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool space( const char c ) {
   return( c == ' ' );
}

/* Split a string on white space returning a vector of strings */

const std::vector<std::string> SplitBetweenBlanks( const std::string& s ) {
   std::vector<std::string> str;
   std::string::const_iterator i = s.begin( );

   while ( i != s.end( ) )
   {
      //look for the next non-blank
      i = std::find_if( i, s.end(), not_space );
      const std::string::const_iterator i2 = std::find_if( i, s.end( ), space );
      str.push_back( std::string( i, i2 ) );
      i = i2;
   }

   return( str );
}

int main(int argc, char ** argv) {

    std::string lat1, lat2;
    std::string arg0;
    std::string arg1;
    std::string arg2;
    std::string line;
    std::vector<std::string> retlines;
    double a1,b1,c1,alpha1,beta1,gamma1,extra1;
    double a2,b2,c2,alpha2,beta2,gamma2,extra2;;
    S6 prim1, prim2;
    std::vector<S6> inputprims;
    double dprim1[6];
    double dprim2[6];
    double rawdist;
    size_t ii,jj,kk;
    char clatsym;
    int argoff;
    int info;

    arg1="";
    arg2="";
    info=0;
    if (argc > 1) { arg1 = std::string(argv[1]); }
    if (argc > 2) { arg2 = std::string(argv[2]); } 
    if (arg1 == "--help" || arg1 == "-h" || arg2 == "--help" || arg2 == "-h") {
        std::cerr
		<< "Usage: cs6dist_mat [--help|-h] print this message and exit" 
		<< std::endl;
        std::cerr
		<< "       cs6dist_mat with no arguments, write distance matrix to cout" 
		<< std::endl;
        std::cerr 
		<< "                                      reading cells from cin" 
		<< std::endl;
        argoff++;
    } 
    if (arg1=="--info" || arg2=="--info")  {
      argoff++;
      info++;
    }


    while (std::getline(std::cin, line)) {
        if (line.size() == 0) break;
        retlines=SplitBetweenBlanks(line);
        if (retlines.size() == 0) break;
        extra1 = 0.;
        lat1 = std::string(retlines[0]);
        clatsym= lat1.substr(0,1)[0];
        if (clatsym == 'D' || clatsym == 'd') {
        extra1 = atof(retlines[7].c_str());
        }
        a1 = atof(retlines[1].c_str());
        b1 = atof(retlines[2].c_str());
        c1 = atof(retlines[3].c_str());
        alpha1 = atof(retlines[4].c_str());
        beta1 = atof(retlines[5].c_str());
        gamma1 = atof(retlines[6].c_str());
        prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1,info);
        inputprims.push_back(prim1);
        if (info) {
          ii = inputprims.size()-1;
          std::cerr << "ii: "<< ii << ": prim1: [" << prim1[0] <<", "<< prim1[1] << ", "<< prim1[2] << ", "
              << prim1[3] << ", " << prim1[4] << ", " << prim1[5] <<  "]" << std::endl;
        }
    }

    {   std::vector<double>  dmat(inputprims.size()*inputprims.size());
        /*std::cerr << "size: " << inputprims.size() << std::endl;*/
        for (ii=0; ii < inputprims.size(); ii++) {
            //std::err << "ii: " << ii << "  " << inputprims[ii] << std::endl;
            prim1 = inputprims[ii];
            for (kk=0; kk < 6; kk++) dprim1[kk] = prim1[kk];
            dmat[ii + ii*inputprims.size()] = 0.;
            for (jj=ii+1; jj < inputprims.size(); jj++) {
                prim2 = inputprims[jj];
                for (kk=0; kk < 6; kk++) dprim2[kk] = prim2[kk];
                /*rawdist = CS6Dist_func(dprim1,dprim2);*/
                rawdist = CS6Dist(dprim1,dprim2);
                dmat[ii+jj*inputprims.size()] = dmat[jj+ii*inputprims.size()] = 0.1*std::sqrt(rawdist);
            }
        }


        for (ii=0; ii < inputprims.size(); ii++) {
            for (jj=0; jj < inputprims.size(); jj++) {
                std::cout <<" "<<dmat[ii+jj*inputprims.size()];
            }
            std::cout << std::endl;
        }
    }

    return 0;
}
