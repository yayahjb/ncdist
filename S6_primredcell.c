
#ifdef __cplusplus

extern "C" {

#endif

#include "S6M_SellingReduce.h"

#include <math.h>
#include <stdlib.h>
#include <ctype.h>


/* Return a primitive reduced S6 cell from cell parameters and a lattice symbol
   testlattice -- 0 (NULL) or a pointer to a C string, the first character if
                  which is a case insensitive lattice type or 'v' for a G6 cell,
                  or d for the first 6 doubles of a D7 cell, or s for an s6 cell.
                    If testlattice is a NULL pointer, p is assumed.
                    If the first character of testlattice is '\0', p is assumed
                    If a lattice type is given the 6 parameters are three cell edge
                      lengths followed by three angles.  if the angles have values
                      between -2*pi and 2*pi, they are assumed to be in radians.
                      Otherwise they are assumed to be in degrees.
   cellparams  -- 0 (NULL) or a pointer to a previously declared double array of
                  at least 6 elements with semantics determined by testlattice
   g6cell      -- 0 (NULL) or a pointer to a previously declared double array of
                  at least 6 elements into which to store the g6 version of the
                  cell defined by testlattice and cellparams
   g6primcell  -- 0 (NULL) or a pointer to a previously declared double array of
                  at least 6 elements into whixh to store the g6 version of the
                  primitive cell, which is not necessarily reduced
   s6primcell  -- 0 (NULL) or a pointer to a previously declared double array of
                  at least 6 elements into whixh to store the s6 version of the
                  primitive cell, which is not necessarily reduced
   s6primred   -- A pointer to a previously declared double array of
                  at least 6 elements into which to store a new S6 primitive
                  S6-reduced cell [b.c, a.c, a.b, a.d, b.d, c.d]
   mat66       -- 0 (NULL) or a pointer to a previously declared double array of
                  at least 36 elements into which to store the matrix elements 
                  (row by row) converting g6 cell into s6primred.
   The return value is a pointer to the values of s6primred if requested and successfully
                  calculated, or 0 (NULL)
 */
double * s6_primredcell(
    char * testlattice,
    double * cellparams,
    double * g6cell, 
    double * g6primcell, 
    double * s6primcell,
    double * s6primred,
    double * mat66) {
    char clatsym;
    double a,b,c,alpha,beta,gamma,extra;
    double xg6cell[6];
    double xcellparams[6]={0.,0.,0.,0.,0.,0.};
    double dmc[36];
    double dg6primcell[6], ds6primcell[6];
    int ii;
    int reduced;
    double * ps6primred;
    if (g6cell == 0) g6cell=xg6cell;
    if (g6primcell == 0) g6primcell=dg6primcell;
    if (s6primcell == 0) s6primcell=ds6primcell;
    if (cellparams == 0) cellparams = xcellparams;

    ps6primred = s6primred;
    if (testlattice==0 || *testlattice=='\0') {
      clatsym='P';
    } else {
      clatsym=*testlattice;
    }
    if (mat66 == 0) mat66 = dmc;
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
            CS6M_CelltoG6(cellparams,g6cell);
            CS6M_LatSymMat66(g6cell,clatsym,mat66,g6primcell);
            break;
        case 'V':
        case 'v':
            for(ii=0; ii<29; ii+=7){
              mat66[ii]=1.;
              mat66[ii+1]=mat66[ii+2]=mat66[ii+3]=mat66[ii+4]=mat66[ii+5]=mat66[ii+6]=0.;
            }
            mat66[35]=1.;
            for(ii=0; ii< 6; ii++) {g6primcell[ii]=cellparams[ii];}
            break;
        case 'D':
        case 'd':
            for(ii=0; ii<29; ii+=7){
              mat66[ii]=1.;
              mat66[ii+1]=mat66[ii+2]=mat66[ii+3]=mat66[ii+4]=mat66[ii+5]=mat66[ii+6]=0.;
            }
            mat66[35]=1.;
            g6primcell[CS6M_G6A2] = cellparams[CS6M_D7A2];
            g6primcell[CS6M_G6B2] = cellparams[CS6M_D7B2];
            g6primcell[CS6M_G6C2] = cellparams[CS6M_D7C2];
            g6primcell[CS6M_G62BC] = cellparams[CS6M_D7_A_D_2]-cellparams[CS6M_D7B2]-cellparams[CS6M_D7C2];
            g6primcell[CS6M_G62AC] = cellparams[CS6M_D7_B_D_2]-cellparams[CS6M_D7A2]-cellparams[CS6M_D7C2];
            g6primcell[CS6M_G62AB] = cellparams[CS6M_D7D2]-cellparams[CS6M_D7A2]-cellparams[CS6M_D7B2]
              -cellparams[CS6M_D7C2]-g6primcell[CS6M_G62BC]-g6primcell[CS6M_G62AC];
           break;
        case 'S':
        case 's':
            for(ii=0; ii<29; ii+=7){
              mat66[ii]=1.;
              mat66[ii+1]=mat66[ii+2]=mat66[ii+3]=mat66[ii+4]=mat66[ii+5]=mat66[ii+6]=0.;
            }
            mat66[35]=1.;
            CS6M_S6toG6(cellparams,g6primcell);
            break;
        default:
            clatsym = 'P';
            CS6M_CelltoG6(cellparams,g6cell);
            CS6M_LatSymMat66(g6cell,clatsym,mat66,g6primcell);
            break;
    }
    CS6M_G6toS6_Mat(g6primcell,s6primcell,mat66);
    reduced=0;
    CS6M_S6Reduce_Mat(s6primcell,s6primred,mat66,reduced);
    if (!reduced) {
      for(ii=0;ii<6;ii++) s6primred[ii]=s6primcell[ii]=0.;
      ps6primred = 0;
    }

    return ps6primred;
}  

#ifdef __cplusplus

}

#endif

