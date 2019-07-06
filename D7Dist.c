//
//  CD7Dist.c
//
//
//  Created by Herbert J. Bernstein and Lawrence C. Andrews on 3/04/17.
//
//
#ifndef D7DIST_C
#define D7DIST_C

#ifdef __cplusplus

extern "C" {

#endif



#include "D7Dist.h"
#include <math.h>

/* Test if outside D7 region
   return 0 if outside */

int d7test(double g[7]) {
    int retval = 0;
    if ( fabs(g[0]+ g[1] + g[2] + g[3] - g[4] - g[5] - g[6]) > 1.e-5 ) retval |= 0x1;
    if ( g[0] + g[3] < g[4] + 1.e-5 ) retval |= 0x2;
    if ( g[1] + g[3] < g[5] + 1.e-5 ) retval |= 0x4;
    if ( g[2] + g[3] < g[6] + 1.e-5 ) retval |= 0x8;
    if ( g[1] + g[2] < g[4] + 1.e-5 ) retval |= 0x10;
    if ( g[0] + g[2] < g[5] + 1.e-5 ) retval |= 0x20;
    if ( g[0] + g[1] < g[6] + 1.e-5 ) retval |= 0x40;
    if ( g[0] <  - 1.e-5 ) retval |= 0x80;
    if ( g[1] <  - 1.e-5 ) retval |= 0x100;
    if ( g[2] <  - 1.e-5 ) retval |= 0x200;
    if ( g[3] <  - 1.e-5 ) retval |= 0x400;
    if ( g[4] <  - 1.e-5 ) retval |= 0x800;
    if ( g[5] <  - 1.e-5 ) retval |= 0x1000;
    if ( g[6] <  - 1.e-5 ) retval |= 0x2000;

    return retval;
}




/*
     Compute the minimal distance between two Delaunay-reduced
     vectors in the D7 Cone following the embedding paths
     to the 7 major boundaries
 */





double D7Dist(double * gvec1,double * gvec2) {
    int rpasses, ir, irt;
    int jr;
    double dist,dist1, dist2, distmin;
    double rgvec1[24][7], rgvec2[24][7];
    double trgvec1[24][7], trgvec2[24][7];
    double ndists[24][24];
    double ndist1[24];
    double ndist2[24];
    double temp[24];
    dist1 = d7minbddist(gvec1);
    dist2 = d7minbddist(gvec2);
    distmin = CD7M_min(dist1,dist2);
    rpasses = NREFL_OUTER_MIN;
    dist = d71234dist(gvec1,gvec2);
    report_double("\n  Entered D7Dist gdist = ",dist,", ");
    report_double_vector("gvec1 = ", gvec1,", ")
    report_double_vector("gvec2 = ", gvec2,";")
    if (dist1+dist2 <  dist*.5 ) {
        rpasses = NREFL_OUTER_MID;
    }
    if (dist1+dist2 <  dist*.1 ) {
        rpasses = NREFL_OUTER_FULL;
    }
    ndists[0][0] = dist = D7Dist_pass(gvec1,gvec2,dist);
/* Collect rpasses-1 transformed vectors */
#pragma omp parallel for schedule(dynamic)
    for (ir = 1; ir < rpasses; ir++) {
        imv7(gvec1,D7Refl[D7Rord[ir]],rgvec1[ir]);
        imv7(gvec2,D7Refl[D7Rord[ir]],rgvec2[ir]);
        ndists[ir][0] = D7Dist_pass(rgvec1[ir],gvec2,dist);
        ndists[0][ir] = D7Dist_pass(gvec1,rgvec2[ir],dist);
    }
    
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (ir = 1; ir < /*rpasses*/ 2; ir++) {
        for (jr = 1; jr < rpasses; jr++) {
            ndists[ir][jr] = D7Dist_pass(rgvec1[ir],rgvec2[jr],dist);
        }
    }
    
#pragma omp flush(dist,ndists)
#pragma omp critical(distminimize)
    for (ir = 0; ir < /*rpasses*/ 2; ir++) {
        for (jr = 0; jr < rpasses; jr++) {
            if (ndists[ir][jr] < dist) dist = ndists[ir][jr];
            report_double_if_changed("\n ndists[ir][jr] ", dist, ", ");
            also_if_changed_report_integer("ir =",ir,", ");
            also_if_changed_report_integer("jr =",jr,"\n");
        }
    }
    return dist;
}



#ifdef __cplusplus

}

#endif
#endif /*D7DIST_C */


