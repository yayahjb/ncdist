#include "G6.h"
#include "Reducer.h"
#include "Delone.h"
#include "Cell.h"
#include "NCDist.h"
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

/*****************************************************************************/

G6 makeprimredprobe( std::string testlattice,
	double a, double b, double c, double alpha, double beta, double gamma )
{
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
    double crootvol;
    Cell rawcell(a,b,c, alpha,beta,gamma);
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
            std::cerr << "Unrecognized lattice symbol "<< testlattice<<" treated as P" << std::endl;
            latsym = "P";
            mc = rawcell.LatSymMat66(latsym);
            primcell = mc*(rawcell.Cell2V6());
            break;
    }
    ret = Reducer::Reduce(primcell,m,redprimcell,0.0);
    ret = Delone::Reduce(primcell,dm,dredprimcell,0.0);
    primredprobe = Cell(redprimcell).CellWithDegrees();
    dprimredprobe = Cell(dredprimcell).CellWithDegrees();
    std::cerr << "Primitive Reduced Probe Cell: " <<
    primredprobe[0]<<" "<<
    primredprobe[1]<<" "<<
    primredprobe[2]<<" "<<
    primredprobe[3]<<" "<<
    primredprobe[4]<<" "<<
    primredprobe[5] << std::endl;
    std::cerr << "Delaunay Primitive Reduced Probe Cell: " <<
    dprimredprobe[0]<<" "<<
    dprimredprobe[1]<<" "<<
    dprimredprobe[2]<<" "<<
    dprimredprobe[3]<<" "<<
    dprimredprobe[4]<<" "<<
    dprimredprobe[5] << std::endl;
    std::cerr << "Volume :" << Cell(redprimcell).Volume() << std::endl;
    crootvol = pow(Cell(redprimcell).Volume(),1./3.);
    Reducer::Reduce((Cell(redprimcell).Inverse()).Cell2V6(),m,reducedBase,0.0);
    recipcell = (Cell(redprimcell).Inverse()).CellWithDegrees();
    
    std::cerr << "Reciprocal of Primitive Probe Cell: " <<
    recipcell[0]<<" "<<
    recipcell[1]<<" "<<
    recipcell[2]<<" "<<
    recipcell[3]<<" "<<
    recipcell[4]<<" "<<
    recipcell[5]<< std::endl;
    std::cerr << "Volume of Reciprocal Cell: " <<
    (Cell(redprimcell).Inverse()).Volume() << std::endl;
    std::cerr << "V7 linearized and scaled: "
    << primredprobe[0]*std::sqrt(6./7.)<<" "
    << primredprobe[1]*std::sqrt(6./7.)<<" "
    << primredprobe[2]*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[0])*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[1])*std::sqrt(6./7.)<<" "
    << std::sqrt(1./reducedBase[2])*std::sqrt(6./7.)<<" "
    << " "<<
    crootvol*std::sqrt(6./7.)<<std::endl;
    if (latsym[0] == 'V' || latsym[0] == 'v') {
        std::cerr << "raw G6 V vector: "
        << primcell[0]<<" "
        << primcell[1]<<" "
        << primcell[2]<<" "
        << primcell[3]<<" "
        << primcell[4]<<" "
        << primcell[5]<<std::endl;
    } else {
        std::cerr << "G6 vector: "
        << primredprobe[0]*primredprobe[0]<<" "
        << primredprobe[1]*primredprobe[1]<<" "
        << primredprobe[2]*primredprobe[2]<<" "
        << 2.*primredprobe[1]*primredprobe[2]*cos(primredprobe[3]*std::atan(1.0)/45.)<<" "
        << 2.*primredprobe[0]*primredprobe[2]*cos(primredprobe[4]*std::atan(1.0)/45.)<<" "
        << 2.*primredprobe[0]*primredprobe[1]*cos(primredprobe[5]*std::atan(1.0)/45.)<<std::endl;
    }
    std::cerr << std::endl;
    return primredprobe;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool not_space( const char c ) {
   return( c != ' ' );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool space( const char c ) {
   return( c == ' ' );
}

/* Split a string on white space returninf a vector of strings */

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
    std::string line;
    std::vector<std::string> retlines;
    double a1,b1,c1,alpha1,beta1,gamma1;
    double a2,b2,c2,alpha2,beta2,gamma2;
    G6 prim1, prim2;
    std::vector<G6> inputprims;
    double rawdist;
    size_t ii,jj,kk;
    char clatsym;
    double dprim1[6];
    double dprim2[6];


    if (argc > 1) arg1 = std::string(argv[1]);
    if (arg1 == "--help" || arg1 == "-h") {
        std::cerr
                << "Usage: ncdist_mat [--help|-h] print this message and exit"
                << std::endl;
        std::cerr
                << "       ncdist_mat with no arguments, write distance matrix to cout"
                << std::endl;
        std::cerr
                << "                                      reading cells from cin"
                << std::endl;
        return -1;
    }

    while (std::getline(std::cin, line)) {
        if (line.size() == 0) break;
        retlines=SplitBetweenBlanks(line);
        if (retlines.size() == 0) break;
        lat1 = std::string(retlines[0]);
        clatsym= lat1.substr(0,1)[0];
        a1 = atof(retlines[1].c_str());
        b1 = atof(retlines[2].c_str());
        c1 = atof(retlines[3].c_str());
        alpha1 = atof(retlines[4].c_str());
        beta1 = atof(retlines[5].c_str());
        gamma1 = atof(retlines[6].c_str());
        prim1 = makeprimredprobe(lat1,a1,b1,c1,alpha1,beta1,gamma1);
        inputprims.push_back(prim1);
        ii = inputprims.size()-1;
        std::cerr << "ii: "<< ii << ": prim1: [" << prim1[0] <<", "<< prim1[1] << ", "<< prim1[2] << ", "
              << prim1[3] << ", " << prim1[4] << ", " << prim1[5] <<  "]" << std::endl;
        std::cerr << "ii: "<< ii << ": inputprims[ii]: [" << inputprims[ii][0] <<", "<<  inputprims[ii][1] << ", "<<  inputprims[ii][2] << ", "
              <<  inputprims[ii][3] << ", " <<  inputprims[ii][4] << ", " <<  inputprims[ii][5] <<  "]" << std::endl;

    }

    { double dmat[inputprims.size()][inputprims.size()];
    std::cout << "size: " << inputprims.size() << std::endl;
    for (ii=0; ii < inputprims.size(); ii++) {
        std::cout << "ii: " << ii << "  " << inputprims[ii] << std::endl;
        prim1 = inputprims[ii];
        for (kk=0; kk < 6; kk++) dprim1[kk] = prim1[kk];
        dmat[ii][ii] = 0.;
        for (jj=ii+1; jj < inputprims.size(); jj++) {
            prim2 = inputprims[jj];
            for (kk=0; kk < 6; kk++) dprim2[kk] = prim2[kk];
            rawdist = NCDist(dprim1,dprim2);
            dmat[ii][jj] = dmat[jj][ii] = 0.1*std::sqrt(rawdist);
        }
    }
    

    for (ii=0; ii < inputprims.size(); ii++) {
        for (jj=0; jj < inputprims.size(); jj++) {
            std::cout <<" "<<dmat[ii][jj];
        }
        std::cout << std::endl;
    }
    }


    return 0;
}

