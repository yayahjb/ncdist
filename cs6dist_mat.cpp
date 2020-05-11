 
#include "G6.h"
#include "D7.h"
#include "S6.h"
#include "Reducer.h"
#include "Delone.h"
#include "LRL_Cell.h"
#include "LRL_Cell_Degrees.h"
#include "D7Dist.h"
/* #include "CS6Dist_func.h" */
#include "CS6Dist.h"
#include "S6M_SellingReduce.h"


#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
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
    std::string arg3;
    std::string arg4;
    std::string line;
    std::vector<std::string> retlines;
    std::vector<std::string> names;
    std::string name;
    int count;
    int nameindex;
    double a1,b1,c1,alpha1,beta1,gamma1,extra1;
    double a2,b2,c2,alpha2,beta2,gamma2,extra2;;
    S6 prim1, prim2;
    G6 g6prim1, g6prim2;
    std::vector<S6> inputprims;
    std::vector<LRL_Cell> inputprimcells;
    LRL_Cell cell1, cell2;
    double dcell1[6], dcell2[6];
    double dprim1[6];
    double dprim2[6];
    double rawdist;
    size_t ii,jj,kk,ll;
    char clatsym;
    int argoff;
    int info;
    int sort;
    int fold;

    arg1="";
    arg2="";
    arg3="";
    arg4="";
    info=0;
    sort=0;
    fold=0;
    if (argc > 1) { arg1 = std::string(argv[1]); }
    if (argc > 2) { arg2 = std::string(argv[2]); } 
    if (argc > 3) { arg3 = std::string(argv[3]); } 
    if (arg1 == "--help" || arg1 == "-h" \
      || arg2 == "--help" || arg2 == "-h" \
      || arg3 == "--help" || arg3 == "-h" \
      || arg4 == "--help" || arg4 == "-h") {
        std::cerr
		<< "Usage: cs6dist_mat [--help|-h] print this message and exit" 
		<< std::endl;
        std::cerr
		<< "       cs6dist_mat [--info] [--sort], write distance matrix to cout" 
		<< std::endl;
        std::cerr 
		<< "                                      reading cells from cin" 
		<< std::endl;
        argoff++;
    } 
    if ( arg1=="--info" || arg2=="--info"|| arg3=="--info" || arg4=="--info" )  {
      argoff++;
      info++;
    }
    if ( arg1=="--sort" || arg2=="--sort"|| arg3=="--sort" || arg4=="--sort" )  {
      argoff++;
      sort++;
    }
    if ( arg1=="--fold" || arg2=="--fold"|| arg3=="--fold" || arg4=="--fold" )  {
      argoff++;
      fold++;
    }


    while (std::getline(std::cin, line)) {
        if (line.size() == 0) break;
        retlines=SplitBetweenBlanks(line);
        if (retlines.size() == 0) break;
        count = retlines.size();
        nameindex = count-1;
        extra1 = 0.;
        lat1 = std::string(retlines[0]);
        clatsym= lat1.substr(0,1)[0];
        name="";
        if (clatsym == 'D' || clatsym == 'd') {
          extra1 = atof(retlines[7].c_str());
          if (nameindex > 7) name = retlines[nameindex];
        } else {
          if (nameindex > 6) name = retlines[nameindex];
        }
        a1 = atof(retlines[1].c_str());
        b1 = atof(retlines[2].c_str());
        c1 = atof(retlines[3].c_str());
        alpha1 = atof(retlines[4].c_str());
        beta1 = atof(retlines[5].c_str());
        gamma1 = atof(retlines[6].c_str());
        prim1 = makeprimredcell(lat1,a1,b1,c1,alpha1,beta1,gamma1,extra1,info);
        inputprims.push_back(prim1);
        names.push_back(name);
        CS6M_S6toG6(prim1,g6prim1);
        CS6M_G6toCell(g6prim1,cell1);
        inputprimcells.push_back(cell1);
        if (info) {
          ii = inputprims.size()-1;
          std::cerr << "ii: "<< ii << ": prim1: [" << prim1[0] <<", "<< prim1[1] << ", "<< prim1[2] << ", "
              << prim1[3] << ", " << prim1[4] << ", " << prim1[5] <<  "]" << " " << name << std::endl;
          std::cerr << "ii: "<< ii << ": g6pim1: [" << g6prim1[0] <<", "<< g6prim1[1] << ", "<< g6prim1[2] << ", "
              << g6prim1[3] << ", " << g6prim1[4] << ", " << g6prim1[5] <<  "]" << " " << name << std::endl;
          std::cerr << "ii: "<< ii << ": cell1: [" << cell1[0] <<", "<< cell1[1] << ", "<< cell1[2] << ", "
              << cell1[3] << ", " << cell1[4] << ", " << cell1[5] <<  "]" << " " << name << std::endl;
        }
    }

    {   std::vector<double>  dmat(inputprims.size()*inputprims.size());
        std::vector<double>  drmsds(inputprims.size()*inputprims.size());
        std::vector<double>  dzs(inputprims.size()*inputprims.size());
        std::vector<size_t>  order(inputprims.size());
        std::vector<double>  sumdist(inputprims.size());
        std::vector<double>  sumdistsq(inputprims.size());
        std::vector<double>  l2norm(inputprims.size());
        std::vector<double>  closecells(inputprims.size());
        int perms[24][4]={
                  {0,1,2,3},
                  {0,1,3,2},
                  {0,2,1,3},
                  {0,2,3,1},
                  {0,3,1,2},
                  {0,3,2,1},
                  {1,0,2,3},
                  {1,0,3,2},
                  {1,2,0,3},
                  {1,2,3,0},
                  {1,3,0,2},
                  {1,3,2,0},
                  {2,1,0,3},
                  {2,1,3,0},
                  {2,0,1,3},
                  {2,0,3,1},
                  {2,3,1,0},
                  {2,3,0,1},
                  {3,1,2,0},
                  {3,1,0,2},
                  {3,2,1,0},
                  {3,2,0,1},
                  {3,0,1,2},
                  {3,0,2,1}};
        /*std::cerr << "size: " << inputprims.size() << std::endl;*/
        for (ii=0; ii < inputprims.size(); ii++) {
            order[ii] = (size_t)ii;
            sumdist[ii] = 0.;
            sumdistsq[ii] = 0.;
            closecells[ii]=1.;
            l2norm[ii] = std::sqrt(inputprims[ii][0]*inputprims[ii][0]
              +inputprims[ii][1]*inputprims[ii][1]
              +inputprims[ii][2]*inputprims[ii][2]
              +inputprims[ii][3]*inputprims[ii][3]
              +inputprims[ii][4]*inputprims[ii][4]
              +inputprims[ii][5]*inputprims[ii][5]);
        }
        for (ii=0; ii < inputprims.size(); ii++) {
            double e3cell1[9], rote3cell2[9];
            //std::err << "ii: " << ii << "  " << inputprims[ii] << std::endl;
            prim1 = inputprims[ii];
            cell1 = inputprimcells[ii];
            for (kk=0; kk < 6; kk++) dprim1[kk] = prim1[kk];
            CS6M_celltoE3(cell1,e3cell1);
            if (info) {
              std::cout << "cell1: ["<<cell1[0]<<","<<cell1[1]<<","<<cell1[2]<<","
              <<cell1[3]<<","<<cell1[4]<<","<<cell1[5]<<"]"<<std::endl;
            }
            dmat[ii + ii*inputprims.size()] = 0.;
            dzs[ii + ii*inputprims.size()] = 0.;;
            for (jj=ii+1; jj < inputprims.size(); jj++) {
                double e3cell2[9];
                double HornMatM[9];
                double HornMatN[16];
                double eigenvec[4];
                double eigenval;
                double rotmat[9];
                double rmsd, rmsdmin;
                double edges[4][3];
                prim2 = inputprims[jj];
                cell2 = inputprimcells[jj];
                for (kk=0; kk < 6; kk++) dprim2[kk] = prim2[kk];
                CS6M_celltoE3(cell2,e3cell2);
                CS6M_E3bycoltoABCD(e3cell2,edges[0],edges[1],edges[2],edges[3]);
                /*rawdist = CS6Dist_func(dprim1,dprim2);*/
                rawdist = CS6Dist(dprim1,dprim2);
                if (rawdist < 1.e-13) rawdist = 0.;
                dmat[ii+jj*inputprims.size()] = dmat[jj+ii*inputprims.size()] = 0.1*std::sqrt(rawdist);
                if (dmat[ii+jj*inputprims.size()] < l2norm[ii]) {
                  closecells[ii]+=1.;
                  sumdist[ii]+=dmat[ii+jj*inputprims.size()];
                  sumdistsq[ii]+=.01*rawdist;
                  closecells[jj]+=1.;
                  sumdist[jj]+=dmat[jj+ii*inputprims.size()];
                  sumdistsq[jj]+=.01*rawdist;
                }
                rmsdmin = 1.e38;
                for (ll=0; ll<24; ll++) {
                  CS6M_vecstoE3bycol(edges[perms[ll][0]],edges[perms[ll][1]],edges[perms[ll][2]],e3cell2) 
                  CS6M_E3cellstoHornMatM(e3cell1,e3cell2,HormMatM);
                  CS6M_HornMatMtoHornMatN(HornMatM,HornMatN);
                  CS6M_EigenMV_Mat44(HornMatN,eigenvec,eigenval);
                  CS6M_QuaterniontoRotmat(eigenvec,rotmat);
                  if (info) {
                    std::cout << "eigenval: " << eigenval << " eigenvec: [" 
                    << eigenvec[0] <<", "<<eigenvec[1] << "," << eigenvec[2] << ", " << eigenvec[3] << "]" <<std::endl;
                  }
                  CS6M_Mat33_Mat33_Mult(rotmat,e3cell2,rote3cell2);
                  CS6M_rmsd(rote3cell2,e3cell1,rmsd);
                  if(rmsdmin > rmsd) { 
                    rmsdmin=rmsd;
                    if (info) {                     
                      std::cout << "perm: "<<perms[ll][0]<<perms[ll][1]<<perms[ll][2]<<perms[ll][3]<< " rmsd: "<<rmsd<<std::endl;
                      std::cout << "HornMatM: ["<<HornMatM[0]<<","<<HornMatM[1]<<","<<HornMatM[2]<<"]"<<std::endl;
                      std::cout << "          ["<<HornMatM[3]<<","<<HornMatM[4]<<","<<HornMatM[5]<<"]"<<std::endl;
                      std::cout << "          ["<<HornMatM[6]<<","<<HornMatM[7]<<","<<HornMatM[8]<<"]"<<std::endl;
                      std::cout << "HornMatN: ["<<HornMatN[0] <<","<<HornMatN[1] <<","<<HornMatN[2] <<","<<HornMatN[3]<<"]"<<std::endl;
                      std::cout << "          ["<<HornMatN[4] <<","<<HornMatN[5] <<","<<HornMatN[6] <<","<<HornMatN[7]<<"]"<<std::endl;
                      std::cout << "          ["<<HornMatN[8] <<","<<HornMatN[9] <<","<<HornMatN[10]<<","<<HornMatN[11]<<"]"<<std::endl;
                      std::cout << "          ["<<HornMatN[12]<<","<<HornMatN[13]<<","<<HornMatN[14]<<HornMatN[15]<<"]"<<std::endl;
                      std::cout << "rotmat: ["<<rotmat[0]<<","<<rotmat[1]<<","<<rotmat[2]<<"]"<<std::endl;
                      std::cout << "        ["<<rotmat[3]<<","<<rotmat[4]<<","<<rotmat[5]<<"]"<<std::endl;
                      std::cout << "        ["<<rotmat[6]<<","<<rotmat[7]<<","<<rotmat[8]<<"]"<<std::endl;

                      std::cout << "rote3cell2 A: [" << rote3cell2[0] << "," << rote3cell2[3] <<  "," << rote3cell2[6] << "]" << std::endl;
                      std::cout << "rote3cell2 B: [" << rote3cell2[1] << "," << rote3cell2[4] <<  "," << rote3cell2[7] << "]" << std::endl;
                      std::cout << "rote3cell2 C: [" << rote3cell2[2] << "," << rote3cell2[5] <<  "," << rote3cell2[8] << "]" << std::endl;
                      std::cout << "rote3cell2 D: [" << -rote3cell2[0]-rote3cell2[1]-rote3cell2[2] << "," 
                                                                         << -rote3cell2[3]-rote3cell2[4]-rote3cell2[5]  <<  "," 
                                                                                                  << rote3cell2[6]-rote3cell2[7]-rote3cell2[8] << "]" << std::endl;
                      std::cout << "   e3cell1 A: [" <<    e3cell1[0] << "," <<    e3cell1[3] <<  "," <<    e3cell1[6] << "]" << std::endl;
                      std::cout << "   e3cell1 B: [" <<    e3cell1[1] << "," <<    e3cell1[4] <<  "," <<    e3cell1[7] << "]" << std::endl;
                      std::cout << "   e3cell1 C: [" <<    e3cell1[2] << "," <<    e3cell1[5] <<  "," <<    e3cell1[8] << "]" << std::endl;
                      std::cout << "   e3cell1 D: [" << -  e3cell1[0]-   e3cell1[1]-   e3cell1[2] << "," 
                                                                         << -   e3cell1[3]-   e3cell1[4]-   e3cell1[5]  <<  "," 
                                                                                                  <<    e3cell1[6]-   e3cell1[7]-   e3cell1[8] << "]" << std::endl;
                      std::cout << "   e3cell2 A: [" <<    e3cell2[0] << "," <<    e3cell2[3] <<  "," <<    e3cell2[6] << "]" << std::endl;
                      std::cout << "   e3cell2 B: [" <<    e3cell2[1] << "," <<    e3cell2[4] <<  "," <<    e3cell2[7] << "]" << std::endl;
                      std::cout << "   e3cell2 C: [" <<    e3cell2[2] << "," <<    e3cell2[5] <<  "," <<    e3cell2[8] << "]" << std::endl;
                      std::cout << "   e3cell2 D: [" << -  e3cell2[0]-   e3cell2[1]-   e3cell2[2] << "," 
                                                                         << -   e3cell2[3]-   e3cell2[4]-   e3cell2[5]  <<  "," 
                                                                                                  <<    e3cell2[6]-   e3cell2[7]-   e3cell2[8] << "]" << std::endl;
                    }
                  }
                }
                if (rmsdmin < 1.e-8) rmsdmin=0.;
                drmsds[ii+jj*inputprims.size()] = drmsds[jj+ii*inputprims.size()] = rmsdmin;
            }
        }
        for (ii=0; ii < inputprims.size(); ii++) {
            for (jj=ii+1; jj < inputprims.size(); jj++) {
                dzs[ii+jj*inputprims.size()] 
                = dmat[ii+jj*inputprims.size()]
                    /std::sqrt((1.e-38
                                +sumdistsq[ii]/(closecells[ii])
                                ));
                dzs[jj+ii*inputprims.size()]
                = dmat[jj+ii*inputprims.size()]
                    /std::sqrt((1.e-38
                                +sumdistsq[jj]/(closecells[jj])
                                ));
            }
        }
        if (sort) {
            size_t gap=inputprims.size()/2;
            int done;
            size_t tmporder;
            done = 0;
            while ( gap > 0 ) {
              done = 1;
              for (jj=0; jj+gap < inputprims.size(); jj+=gap) {
                  if ( sumdistsq[order[jj]]*closecells[order[jj+gap]] > sumdistsq[order[jj+gap]]*closecells[order[jj]] ) {
                      done = 0;
                      tmporder = order[jj];
                      order[jj] = order[jj+gap];
                      order[jj+gap] = tmporder; 
                  }
              }
              if (done) gap /= 2;
            }
        }

        std::cout << "cs6dist: ";
        for (jj=0; jj < names.size(); jj++) {
          std::cout << " " << std::setw(12) << names[order[jj]];
          if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
        }
        std::cout << std::endl;
        for (ii=0; ii < inputprims.size(); ii++) {
            std::cout <<  " " << std::setw(12) <<names[order[ii]];
            for (jj=0; jj < inputprims.size(); jj++) {
                std::cout <<" "<< std::setw(12) << dmat[order[ii]+order[jj]*inputprims.size()];
                if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
            }
            std::cout << std::endl;
        }

        std::cout << "z-scores: ";
        for (jj=0; jj < names.size(); jj++) {
          std::cout << " " << std::setw(12) << names[order[jj]];
          if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
        }
        std::cout << std::endl;
        for (ii=0; ii < inputprims.size(); ii++) {
            std::cout <<  " " << std::setw(12) <<names[order[ii]];
            for (jj=0; jj < inputprims.size(); jj++) {
                std::cout <<" "<< std::setw(12) << dzs[order[ii]+order[jj]*inputprims.size()];
                if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
            }
            std::cout << std::endl;
        }

        std::cout << "   rmsd: ";
        for (jj=0; jj < names.size(); jj++) {
          std::cout << " " << std::setw(12) << names[order[jj]];
          if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
        }
        std::cout << std::endl;
        for (ii=0; ii < inputprims.size(); ii++) {
            std::cout <<  " " << std::setw(12) << names[order[ii]];
            for (jj=0; jj < inputprims.size(); jj++) {
                std::cout <<" "<< std::setw(12) << drmsds[order[ii]+order[jj]*inputprims.size()];
                if ( fold && (jj+1)%8 ==0 ) {std::cout << std::endl;}
            }
            std::cout << std::endl;
        }
    }

    return 0;
}
