//
//  CD7Dist.c
//
//
//  Created by Herbert J. Bernstein and Lawrence C. Andrews on 3/04/17.
//
//

/* The projectors for the 15 base types (5-D boundaries
 in G6) plus a few extra for internal boundaries
 Note that the array indices are swapped from the
 Fortan versions */
#ifndef D7DIST_H
#define D7DIST_H


/* #define NCDIST_DEBUG */
#define NCDIST_NO_OUTER_PASS

#include <math.h>
#include <float.h>

#ifdef NCDIST_DEBUG
static double oldvalue;
#include <cstdio>
#define d7report_double(prolog, value, epilog) \
oldvalue=value; fprintf(stderr, "%s%g%s", prolog, value, epilog);
#define d7report_integer(prolog, value, epilog) \
fprintf(stderr, "%s%d%s", prolog, value, epilog);
#define d7report_double_if_changed(prolog, value, epilog) \
changed=0; if (fabs(value-oldvalue)>1.e-8*(fabs(value)+fabs(oldvalue)+1.e-12)) {oldvalue=value; changed=1; fprintf(stderr, "%s%g%s", prolog, value, epilog);}
#define also_if_changed_report(prolog, value, epilog) \
if(changed) {fprintf(stderr, "%s%s%s", prolog, value, epilog);}
#define also_if_changed_d7report_integer(prolog, value, epilog) \
if(changed) {fprintf(stderr, "%s%d%s", prolog, value, epilog);}
#define also_if_changed_d7report_double(prolog, value, epilog) \
if(changed) {fprintf(stderr, "%s%g%s", prolog, value, epilog);}
#define also_if_changed_d7report_double_vector(prolog, value, epilog) \
if(changed) {fprintf(stderr, "%s[%g  %g  %g  %g  %g  %g]%s", prolog, value[0], value[1], value[2], value[3], value[4], value[5], epilog);}
#define d7report_double_vector(prolog, value, epilog) \
{fprintf(stderr, "%s[%g  %g  %g  %g  %g  %g]%s", prolog, value[0], value[1], value[2], value[3], value[4], value[5], epilog);}
#else
#define d7report_double(prolog, value, epilog)
#define d7report_integer(prolog, value, epilog)
#define d7report_double_if_changed(prolog, value, epilog)
#define also_if_changed_report(prolog, value, epilog)
#define also_if_changed_d7report_integer(prolog, value, epilog)
#define also_if_changed_d7report_double(prolog, value, epilog)
#define also_if_changed_d7report_double_vector(prolog, value, epilog)
#define d7report_double_vector(prolog, value, epilog)
#endif





#define CD7M_min(a,b) (((a)<(b))?(a):(b))
#define CD7M_max(a,b) (((a)<(b))?(b):(a))

#define ND7BND 9

#define D7P_1 0
#define D7P_2 1
#define D7P_3 2
#define D7P_4 3
#define D7P_5 4
#define D7P_6 5
#define D7P_7 6
#define D7P_8 7
#define D7P_9 8

#define ND7BDPRJ 9



#define P_1   0
#define P_2   1
#define P_3   2
#define P_4   3
#define P_5   4
#define P_6   5
#define P_7   6
#define P_8   7
#define P_9   8
#define P_A   9
#define P_B  10
#define P_C  11
#define P_D  12
#define P_E  13
#define P_F  14
#define P_6C 15
#define P_67 16
#define P_9A 17
#define P_CD 18
#define P_12 19
#define P_8B 20
#define P_8E 21
#define P_8F 22
#define P_BF 23
#define P_EF 24
#define P_28F 25
#define P_2BF 26
#define P_2EF 27
#define P_269 28
#define P_26C 29
#define P_2F  30
#define P_27  31
#define P_2A  32
#define P_2D  33
#define P_28E  34
#define P_28B  35


/* D7Dist Boundary projectors */

static double d7prj[ND7BDPRJ][49]= {
    /*prj[D7P_1]  ,1,*/
    {0.3571428571428572,0.3571428571428572,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        0.3571428571428572,0.3571428571428572,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,0.8571428571428571,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        0.8571428571428571,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,0.8571428571428571,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,0.8571428571428571,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        0.8571428571428571},
    /*prj[D7P_2]   2,*/
    {0.8571428571428571,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,0.3571428571428572,0.3571428571428572,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,0.3571428571428572,0.357142571428572,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        0.8571428571428571,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,0.8571428571428571,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,0.8571428571428571,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        0.8571428571428571},
    /*prj[D7P_3]   3,*/
    {0.8571428571428571,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,0.8571428571428571,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,0.3571428571428572,
        0.3571428571428572,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,0.3571428571428572,
        0.3571428571428572,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,0.8571428571428571,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,0.8571428571428571,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        0.8571428571428571},
    /*prj[D7P_4]   4,*/
    {0.75,0.0,0.0,-0.25,0.0,0.25,0.25,
        0.0,0.6666666666666666,-0.3333333333333333,0.0,
        0.3333333333333333,0.0,0.0,
        0.0,-0.3333333333333333,0.6666666666666666,0.0,
        0.3333333333333333,0.0,0.0,-0.25,0.0,0.0,0.75,0.0,0.25,0.25,
        0.0,0.3333333333333333,0.3333333333333333,0.0,
        0.6666666666666666,0.0,0.0,0.25,0.0,0.0,0.25,0.0,0.75,-0.25,
        0.25,0.0,0.0,0.25,0.0,-0.25,0.75},
    /*prj[D7P_5]   5 */
    {0.6666666666666666,0.0,0.0,-0.3333333333333333,
        0.3333333333333333,0.0,0.0,0.0,0.75,-0.25,0.0,0.0,0.25,0.25,
        0.0,-0.25,0.75,0.0,0.0,0.25,0.25,
        -0.3333333333333333,0.0,0.0,0.6666666666666666,
        0.3333333333333333,0.0,0.0,
        0.3333333333333333,0.0,0.0,0.3333333333333333,
        0.6666666666666666,0.0,0.0,0.0,0.25,0.25,0.0,0.0,0.75,-0.25,
        0.0,0.25,0.25,0.0,0.0,-0.25,0.75},
    /*prj[D7P_6]   6 */
    {0.6666666666666666,0.0,-0.3333333333333333,0.0,0.0,
        0.3333333333333333,0.0,0.0,0.75,0.0,-0.25,0.25,0.0,0.25,
        -0.3333333333333333,0.0,0.6666666666666666,0.0,0.0,
        0.3333333333333333,0.0,0.0,-0.25,0.0,0.75,0.25,0.0,0.25,
        0.0,0.25,0.0,0.25,0.75,0.0,-0.25,
        0.3333333333333333,0.0,0.3333333333333333,0.0,0.0,
        0.6666666666666666,0.0,0.0,0.25,0.0,0.25,-0.25,0.0,0.75},
    /*prj[D7P_7]   7 */
    {0.75,0.0,-0.25,0.0,0.25,0.0,0.25,
        0.0,0.6666666666666666,0.0,-0.3333333333333333,0.0,
        0.3333333333333333,0.0,-0.25,0.0,0.75,0.0,0.25,0.0,0.25,
        0.0,-0.3333333333333333,0.0,0.6666666666666666,0.0,
        0.3333333333333333,0.0,0.25,0.0,0.25,0.0,0.75,0.0,-0.25,
        0.0,0.3333333333333333,0.0,0.3333333333333333,0.0,
        0.6666666666666666,0.0,0.25,0.0,0.25,0.0,-0.25,0.0,0.75},
    /*prj[D7P_8]   8 */
    {0.6666666666666666,-0.3333333333333333,0.0,0.0,0.0,0.0,
        0.3333333333333333,
        -0.3333333333333333,0.6666666666666666,0.0,0.0,0.0,0.0,
        0.3333333333333333,0.0,0.0,0.75,-0.25,0.25,0.25,0.0,
        0.0,0.0,-0.25,0.75,0.25,0.25,0.0,
        0.0,0.0,0.25,0.25,0.75,-0.25,0.0,
        0.0,0.0,0.25,0.25,-0.25,0.75,0.0,
        0.3333333333333333,0.3333333333333333,0.0,0.0,0.0,0.0,
        0.6666666666666666},
    /*prj[D7P_9]   9 */
    {0.75,-0.25,0.0,0.0,0.25,0.25,0.0,
        -0.25,0.75,0.0,0.0,0.25,0.25,0.0,
        0.0,0.0,0.6666666666666666,-0.3333333333333333,0.0,0.0,
        0.3333333333333333,
        0.0,0.0,-0.3333333333333333,0.6666666666666666,0.0,0.0,
        0.3333333333333333,0.25,0.25,0.0,0.0,0.75,-0.25,0.0,
        0.25,0.25,0.0,0.0,-0.25,0.75,0.0,
        0.0,0.0,0.3333333333333333,0.3333333333333333,0.0,0.0,
        0.6666666666666666}};

static double d7prj_perp[ND7BDPRJ][49]= {
    /*prj_perp[D7P_1]   1 */
    {0.6428571428571429,-0.3571428571428572,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        -0.3571428571428572,0.6428571428571429,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428},
    /*prj_perp[D7P_2]   2 */
    {0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.6428571428571429,-0.3571428571428572,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,-0.3571428571428572,0.6428571428571429,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428},
    /*prj_perp[D7P_3]   3 */
    {0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,0.6428571428571429,
        -0.3571428571428572,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        0.1428571428571428,0.1428571428571428,-0.3571428571428572,
        0.6428571428571429,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428,
        -0.1428571428571428,-0.1428571428571428,-0.1428571428571428,
        -0.1428571428571428,0.1428571428571428,0.1428571428571428,
        0.1428571428571428},
    /*prj_perp[D7P_4]   4 */
    {0.25,0.0,0.0,0.25,0.0,-0.25,-0.25,
        0.0,0.3333333333333333,0.3333333333333333,0.0,
        -0.3333333333333333,0.0,0.0,
        0.0,0.3333333333333333,0.3333333333333333,0.0,
        -0.3333333333333333,0.0,0.0,
        0.25,0.0,0.0,0.25,0.0,-0.25,-0.25,
        0.0,-0.3333333333333333,-0.3333333333333333,0.0,
        0.3333333333333333,0.0,0.0,
        -0.25,0.0,0.0,-0.25,0.0,0.25,0.25,
        -0.25,0.0,0.0,-0.25,0.0,0.25,0.25},
    /*prj_perp[D7P_5]   5 */
    {0.3333333333333333,0.0,0.0,0.3333333333333333,
        -0.3333333333333333,0.0,0.0,
        0.0,0.25,0.25,0.0,0.0,-0.25,-0.25,
        0.0,0.25,0.25,0.0,0.0,-0.25,-0.25,
        0.3333333333333333,0.0,0.0,0.3333333333333333,
        -0.3333333333333333,0.0,0.0,
        -0.3333333333333333,0.0,0.0,-0.3333333333333333,
        0.3333333333333333,0.0,0.0,
        0.0,-0.25,-0.25,0.0,0.0,0.25,0.25,
        0.0,-0.25,-0.25,0.0,0.0,0.25,0.25},
    /*prj_perp[D7P_6]   6 */
    {0.3333333333333333,0.0,0.3333333333333333,0.0,0.0,
        -0.3333333333333333,0.0,0.0,0.25,0.0,0.25,-0.25,0.0,-0.25,
        0.3333333333333333,0.0,0.3333333333333333,0.0,0.0,
        -0.3333333333333333,0.0,0.0,0.25,0.0,0.25,-0.25,0.0,-0.25,
        0.0,-0.25,0.0,-0.25,0.25,0.0,0.25,
        -0.3333333333333333,0.0,-0.3333333333333333,0.0,0.0,
        0.3333333333333333,0.0,0.0,-0.25,0.0,-0.25,0.25,0.0,0.25},
    /*prj_perp[D7P_7]   7 */
    {0.25,0.0,0.25,0.0,-0.25,0.0,-0.25,
        0.0,0.3333333333333333,0.0,0.3333333333333333,0.0,
        -0.3333333333333333,0.0,0.25,0.0,0.25,0.0,-0.25,0.0,-0.25,
        0.0,0.3333333333333333,0.0,0.3333333333333333,0.0,
        -0.3333333333333333,0.0,-0.25,0.0,-0.25,0.0,0.25,0.0,0.25,
        0.0,-0.3333333333333333,0.0,-0.3333333333333333,0.0,
        0.3333333333333333,0.0,-0.25,0.0,-0.25,0.0,0.25,0.0,0.25},
    /*prj_perp[D7P_8]   8 */
    {0.3333333333333333,0.3333333333333333,0.0,0.0,0.0,0.0,
        -0.3333333333333333,
        0.3333333333333333,0.3333333333333333,0.0,0.0,0.0,0.0,
        -0.3333333333333333,0.0,0.0,0.25,0.25,-0.25,-0.25,0.0,
        0.0,0.0,0.25,0.25,-0.25,-0.25,0.0,
        0.0,0.0,-0.25,-0.25,0.25,0.25,0.0,
        0.0,0.0,-0.25,-0.25,0.25,0.25,0.0,
        -0.3333333333333333,-0.3333333333333333,0.0,0.0,0.0,0.0,
        0.3333333333333333},
    /*prj_perp[D7P_9]   9 */
    {0.25,0.25,0.0,0.0,-0.25,-0.25,0.0,
        0.25,0.25,0.0,0.0,-0.25,-0.25,0.0,
        0.0,0.0,0.3333333333333333,0.3333333333333333,0.0,0.0,
        -0.3333333333333333,
        0.0,0.0,0.3333333333333333,0.3333333333333333,0.0,0.0,
        -0.3333333333333333,-0.25,-0.25,0.0,0.0,0.25,0.25,0.0,
        -0.25,-0.25,0.0,0.0,0.25,0.25,0.0,
        0.0,0.0,-0.3333333333333333,-0.3333333333333333,0.0,0.0,
        0.3333333333333333}};



/* The following matrices are the transformation
 matrices that may be applied at the associated
 boundaries.  *** An addition post permutation may
 be needed to preserve D7 ordering *** */

/* D7 post permuations {1,2,3,4,5,6,7},   Ident
                       {1,2,4,3,6,5,7},   D7M_3
                       {1,3,2,4,5,7,6},   D7M_2
                       {1,3,4,2,7,5,6},   D7M_1.D7M_2
                       {1,4,2,3,6,7,5},
                       {1,4,3,2,7,6,5},
                       {3,1,2,4,7,5,6},
                       {4,1,2,3,7,6,5}
 */

#define D7Refl_1 0
#define D7Refl_2 1
#define D7Refl_3 2
#define D7Refl_4 3
#define D7Refl_term -1

static int D7Refl[4][49] = {
    /* D7Refl_1 {1,2,3,4,5,6,7},   Ident */
    {1,0,0,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,0,0,1},
    /* D7Refl_2 {2,1,3,4,6,5,7},   a <-> b */
    {0,1,0,0,0,0,0,
        1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1},
    /* D7Refl_3 {1,3,2,4,5,7,6},   b <-> c */
    {1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1,
        0,0,0,0,0,1,0},
    /* D7Refl_4 {1,2,4,3,6,5,7},   c <-> d */
    {1,0,0,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1},
    
};

/*  The 24 D7 permuations, the first 8 of which are seen in the
    boundary transforms */

static int D7Perm[24][8] = {
    
    /* D7Perm_1 {1,2,3,4,5,6,7},   Ident   OK */
    {D7Refl_1,D7Refl_term,0,0,0,0,0,0},
    /* D7Perm_2 {1,2,4,3,6,5,7}            OK */
    {D7Refl_4,D7Refl_term,0,0,0,0,0,0},
    /* D7Perm_3 {1,3,2,4,5,7,6}            OK */
    {D7Refl_3,D7Refl_term,0,0,0,0,0,0},
    /* D7Perm_4 {1,3,4,2,7,5,6},           OK */
    {D7Refl_4,D7Refl_3,D7Refl_term,0,0,0,0,0},
    /* D7Perm_5 {1,4,2,3,6,7,5}            OK */
    {D7Refl_3,D7Refl_4,D7Refl_term,0,0,0,0,0},
    /* D7Perm_6  {1,4,3,2,7,6,5}           OK */
    {D7Refl_4,D7Refl_3,D7Refl_4,D7Refl_term,0,0,0,0},
    /* D7Perm_7  {3,1,2,4,7,5,6}           OK */
    {D7Refl_2,D7Refl_3,D7Refl_term,0,0,0,0,0},
    /* D7Perm_8  {4,1,2,3,7,6,5}           OK */
    {D7Refl_2,D7Refl_3,D7Refl_4,D7Refl_term,0,0,0,0},
    /* D7Perm_9  {2, 1, 3, 4}              OK */
    {D7Refl_2,D7Refl_term,0,0,0,0,0,0},
    /* D7Perm_10 {2, 1 ,4, 3}              OK */
    {D7Refl_4,D7Refl_2,D7Refl_term,0,0,0,0,0},
    /* D7Perm_11 {2, 3, 1, 4}              OK */
    {D7Refl_3,D7Refl_2,D7Refl_term,0,0,0,0,0},
    /* D7Perm_12 {2, 3, 4, 1}              OK */
    {D7Refl_4,D7Refl_3,D7Refl_2,D7Refl_term,0,0,0,0},
    /* D7Perm_13 {2, 4, 1, 3}              OK */
    {D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_term,0,0,0,0},
    /* D7Perm_14 {2, 4, 3, 1}              OK */
    {D7Refl_4,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_term,0,0,0},
    /* D7Perm_15 {3, 1, 4, 2}              OK */
    {D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_term,0,0,0,0},
    /* D7Perm_16 {3, 2, 1, 4}              OK */
    {D7Refl_3,D7Refl_2,D7Refl_3,D7Refl_term,0,0,0,0},
    /* D7Perm_17 {3, 2, 4, 1}              OK */
    {D7Refl_4,D7Refl_3,D7Refl_2,D7Refl_3,D7Refl_term,0,0,0},
    /* D7Perm_18 {3, 4, 1, 2}              OK */
    {D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_term,0,0,0},
    /* D7Perm_19 {3, 4, 2, 1}              OK */
    {D7Refl_4,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_term,0,0},
    /* D7Perm_20 {4, 1, 3, 2}              OK */
    {D7Refl_2,D7Refl_4,D7Refl_3,D7Refl_4,D7Refl_term,0,0,0},
    /* D7Perm_21 {4, 2, 1, 3}              OK */
    {D7Refl_2,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_term,0,0,0},
    /* D7Perm_22 {4, 2, 3, 1}              OK */
    {D7Refl_2,D7Refl_4,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_term,0,0},
    /* D7Perm_23 {4, 3, 1, 2}              OK */
    {D7Refl_2,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_term,0,0},
    /* D7Perm_24 {4, 3, 2, 1}              OK */
    {D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_4,D7Refl_2,D7Refl_3,D7Refl_term,0}
};

static int D7MS[9][49] = {
    
    /* D7M_1 d1 <-> d2, d5 <-> d6 {2,1,3,4,6,5,7} */
    {0,1,0,0,0,0,0,
        1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1},

    /* D7M_2 d2 <-> d3, d6 <-> d7 {1,3,2,4,5,7,6} */
    {1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1,
        0,0,0,0,0,1,0},

    /* D7M_3 d3 <-> d4, d5 <-> d6 {1,2,4,3,6,5,7}*/
    {1,0,0,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1},
    
    /* D7M_4 d1->d2, d2->d3, d3->d6, d4->d7, d5->d1, d6->d4, d7->2d2+2d3-d5
       and post permutations 1,2,3,4,5,6,7,8
     */
    {0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,0,0,1,
        1,0,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,2,2,0,-1,0,0},
    
    /* D7M_5 d2->d4, d3->d5, d4->d1, d5->d2, d6->d3, d7->-2d2-2d3+d5+2d6+2d7
       and post permutations 1,2,3,4,5,6,7,8 */
    {1,0,0,0,0,0,0,
        0,0,0,1,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,0,0,1,
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,-2.-2,0,1,2,2},
    
    /* D7M_6 d2->d3, d3->d5, d4->d1, d5->d2, d6->d4, d6-2d1+2d3-d6
       and post permutations 1,2,3,4,5,6,7,8 */
    
    {1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        2,0,2,0,0,-1,0},

    /* D7M_7 d2->d3, d3->d5, d4->d1, d5->d2, d6->d4, d6-2d1+2d3-d6
     and post permutations 1,2,3,4,5,6,7,8, all are case 6 cells */
    
    {1,0,0,0,0,0,0,
        0,0,1,0,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,0,1,
        0,1,0,0,0,0,0,
        0,0,0,1,0,0,0,
        2,0,2,0,0,-1,0},
    
    /* D7M_8 d3->d5, d4->d6, d5->d3, d6->d4, d7->2d1+2d2-d7
     and post permutations 1,2,3,5,8,7 */
    
    {1,0,0,0,0,0,0,
        0,1,0,0,0,0,0,
        0,0,0,0,1,0,0,
        0,0,0,0,0,1,0,
        0,0,1,0,0,0,0,
        0,0,0,1,0,0,0,
        2,2,0,0,0,0,-1},
    
    /* D7M_9 d1->d3, d2->d6, d3->d5, d5->-2d1-2d2+2d5+2d6+d7,d6->d2,d7->d1 */
    
    {0,0,1,0,0,0,0,
        0,0,0,0,0,1,0,
        0,0,0,0,1,0,0,
        0,0,0,1,0,0,0,
        -2,-2,0,0,2,2,1,
        0,1,0,0,0,0,0,
        1,0,0,0,0,0,0}

};




#undef M_E
#define M_1 0
#define M_2 1
#define M_3 2
#define M_4 3
#define M_5 4
#define M_6 5
#define M_7 6
#define M_8 7
#define M_9 8
#define M_A 9
#define M_B 10
#define M_C 11
#define M_D 12
#define M_E 13
#define M_F 14

#define bdf_P_2F 0
#define bdf_M_2_P_2F  1           /* on 2F */
#define bdf_M_F_M_2_P_2F 2        /* on 2F */
#define bdf_M_F_P_2F 3            /* on 2F */
#define bdf_M_2_M_F_P_2F 4        /* on 2F */

#define bdf_P_69 1                /* on 69 */
#define bdf_M_6_P_69 5            /* on 8F (+++--- tunnel) */
#define bdf_M_9_P_69 6            /* on BF (+++--- tunnel) */
#define bdf_M_F_M_6_P_69 7        /* in BF (+++--- tunnel) */
#define bdf_M_F_M_9_P_69 8        /* in 8F (+++--- tunnel) */

#define bdf_P_6C 2                /* on 6C */
#define bdf_M_C_P_6C 9            /* on EF (+++--- tunnel) */
#define bdf_M_F_M_C_P_6C 10       /* on EF (+++--- tunnel) */
#define bdf_M_6_P_6C 11           /* on 8E (+++--- tunnel) */

#define bdf_P_8B  3

#define bdf_P_8E  4
#define bdf_M_8_P_8E 12           /* on 8B */

#define bdf_P_8F 5                /* on 8F */
#define bdf_M_8_P_8F  13          /* on 69 (---+++ tunnel) */
#define bdf_M_F_P_8F 14           /* on BF */
#define bdf_M_B_M_F_P_8F  15      /* on 69 (---+++ tunnel) */

#define bdf_P_BF 6                /* on BF */
#define bdf_M_B_P_BF 16           /* on 69 (---+++ tunnel) */
#define bdf_M_F_P_BF 17           /* on 8F */
#define bdf_M_8_M_F_P_BF 18       /* on 69 (---+++ tunnel)*/

#define bdf_P_EF 7                /* on EF */
#define bdf_M_E_P_EF 19           /* on 6C (---+++ tunnel) */
#define bdf_M_6_M_E_P_EF 20       /* on 8E */
#define bdf_M_F_P_EF 21           /* on EF */

#define bdf_P_269 8
#define bdf_M_2_P_269 22          /* on 26C */
#define bdf_M_6_P_269 23          /* on 28F */
#define bdf_M_9_P_269 24          /* on 2BF */
#define bdf_M_2_M_6_P_269 25      /* on 28F */
#define bdf_M_2_M_9_P_269 26      /* on 2EF */


#define bdf_P_26C 9
#define bdf_M_2_P_26C 27          /* on 269 */
#define bdf_M_6_P_26C 28          /* on 28E */
#define bdf_M_C_P_26C 29          /* on 2EF */
#define bdf_M_2_M_6_P_26C 30      /* on 28B */
#define bdf_M_2_M_C_P_26C 31      /* on 2BF */

#define bdf_P_28F 10
#define bdf_M_2_P_28F 32          /* on 28F */
#define bdf_M_8_P_28F 33          /* on 269 */
#define bdf_M_F_P_28F 34          /* on 2BF */
#define bdf_M_F_M_2_P_28F 35      /* on 2BF */
#define bdf_M_2_M_F_M_2_P_28F 36  /* on 2EF */



#define bdf_P_2BF 11
#define bdf_M_2_P_2BF  37         /* on 2EF */
#define bdf_M_B_P_2BF  38         /* on 28F */
#define bdf_M_F_P_2BF  39         /* on 2BF */
#define bdf_M_F_M_2_P_2BF 40      /* on 2EF */
#define bdf_M_2_M_F_M_2_P_2BF 41  /* on 2BF */



#define bdf_P_2EF 12
#define bdf_M_2_P_2EF 42          /* on 2BF */
#define bdf_M_E_P_2EF 43          /* on 26C */
#define bdf_M_F_P_2EF 44          /* on 2EF */
#define bdf_M_F_M_2_P_2EF 45      /* on 28F */
#define bdf_M_2_M_F_M_2_P_2EF 46  /* on 28F */

#define bdf_P_27 13
#define bdf_M_2_P_27  47          /* on 27 */
#define bdf_M_7_M_2_P_27 48       /* on 27 */
#define bdf_M_7_P_27 49           /* on 27 */
#define bdf_M_2_M_7_P_27 50       /* on 27 */

#define bdf_P_2A 14
#define bdf_M_2_P_2A 51           /* on 2D */
#define bdf_M_D_M_2_P_2A 52       /* on 2D */
#define bdf_M_A_P_2A 53           /* on 2A */
#define bdf_M_2_M_A_P_2A 54       /* on 2D */

#define bdf_P_2D 15
#define bdf_M_2_P_2D 55           /* on 2A */
#define bdf_M_A_M_2_P_2D 56       /* on 2A */
#define bdf_M_D_P_2D 57           /* on 2D */
#define bdf_M_2_M_D_P_2D 58       /* on 2A */

#define bdf_P_28E 16
#define bdf_M_2_P_28E  59         /* on 28B */
#define bdf_M_8_P_28E  60         /* on 26C */
#define bdf_M_2_M_8_P_28E 61      /* on 269 */

#define bdf_P_28B 17
#define bdf_M_2_P_28B 62          /* on 28E */
#define bdf_M_8_M_2_P_28B 63      /* on 26C */
#define bdf_M_2_M_8_M_2_P_28B 64  /* on 269 */


#define NMPGS            65
#define NPGS             18




/* Compute the dot product of 2 D7 vectors */

static double d7dotprod(double v1[7], double v2[7]) {
    
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3]+v1[4]*v2[4]+v1[5]*v2[5]+v1[6]*v2[6];
    
}


/* Compute the simple Euclidean distance between 2 D7 vectors */

static double d7eucldistsq(double v1[7], double v2[7]) {
    
    double vtemp;
    int ii;
    double distsq;
    
    distsq = 0.;
    
    for (ii = 0; ii < 7; ii++ ) {
        vtemp = v1[ii]-v2[ii];
        distsq += vtemp*vtemp;
    }
    return distsq;
}

static double d7eucldist(double v1[7], double v2[6]) {
    
    return sqrt(d7eucldistsq(v1,v2));
    
}



#define CD7M_min2(a,b) (((a)<(b))?(a):(b))
#define CD7M_min3(a,b,c) ((CD7M_min2(a,b)<(c))?(CD7M_min2(a,b)):(c))
#define CD7M_min4(a,b,c,d) ((CD7M_min3(a,b,c)<(d))?(CD7M_min3(a,b,c)):(d))
#define CD7M_min5(a,b,c,d,e) ((CD7M_min4(a,b,c,d)<(e))?(CD7M_min4(a,b,c,d)):(e))
#define CD7M_min6(a,b,c,d,e,f) ((CD7M_min5(a,b,c,d,e)<(f))?(CD7M_min5(a,b,c,d,e)):(f))
#define CD7M_min7(a,b,c,d,e,f,g) ((CD7M_min6(a,b,c,d,e,f)<(g))?(CD7M_min6(a,b,c,d,e,f)):(g))
#define CD7M_min8(a,b,c,d,e,f,g,h) ((CD7M_min7(a,b,c,d,e,f,g)<(h))?(CD7M_min7(a,b,c,d,e,f,g)):(h))


/* Macro versions of distances */

#define CD7M_normsq(v) v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]+v[4]*v[4]+v[5]*v[5]
#define CD7M_norm(v) sqrt(CD7M_normsq(v))
#define CD7M_d7normsq(v) v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]+v[4]*v[4]+v[5]*v[5]+v[6]*v[6]
#define CD7M_d7norm(v) sqrt(CD7M_d7normsq(v))

#define CD7M_gdistsq(v1,v2) \
( \
(v1[0]-v2[0])*(v1[0]-v2[0])+\
(v1[1]-v2[1])*(v1[1]-v2[1])+\
(v1[2]-v2[2])*(v1[2]-v2[2])+\
(v1[3]-v2[3])*(v1[3]-v2[3])+\
(v1[4]-v2[4])*(v1[4]-v2[4])+\
(v1[5]-v2[5])*(v1[5]-v2[5]))

#define CD7M_gdist(v1,v2) sqrt(CD7M_gdistsq(v1,v2))

#define CD7M_d7eucldistsq(v1,v2) \
( \
(v1[0]-v2[0])*(v1[0]-v2[0])+\
(v1[1]-v2[1])*(v1[1]-v2[1])+\
(v1[2]-v2[2])*(v1[2]-v2[2])+\
(v1[3]-v2[3])*(v1[3]-v2[3])+\
(v1[4]-v2[4])*(v1[4]-v2[4])+\
(v1[5]-v2[5])*(v1[5]-v2[5])+\
(v1[6]-v2[6])*(v1[6]-v2[6])\
)

#define CD7M_d7eucldist(v1,v2) sqrt(CD7M_d7eucldistsq(v1,v2))


/*   Macro version of g456dist
 Compute the best distance between 2 G6 vectors
 allowing for cell-preserving sign changes in
 g4,5,6
 */

#define CD7M_g456distsq(v1,v2) \
fabs( CD7M_gdistsq(v1,v2)+\
 4.*CD7M_min(CD7M_min(CD7M_min(0.,       \
                            v1[3]*v2[3]+v1[4]*v2[4]), \
                    v1[3]*v2[3]+v1[5]*v2[5]), \
             v1[4]*v2[4]+v1[5]*v2[5]))

#define CD7M_g456dist(v1,v2) sqrt(CD7M_g456distsq(v1,v2))

/*     Compute the best distance between 2 G6 vectors
 allowing for permulations of g1, g2, g3 as
 well as sign changes
 */

#define CD7M_g456distsq_byelem(v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26) \
fabs((v11-v21)*(v11-v21)+(v12-v22)*(v12-v22)+(v13-v23)*(v13-v23) + \
(v14-v24)*(v14-v24)+(v15-v25)*(v15-v25)+(v16-v26)*(v16-v26) + \
4.*CD7M_min(CD7M_min(CD7M_min(0.,v14*v24+v15*v25),v14*v24+v16*v26),v15*v25+v16*v26))



/*     Compute the best distance between 2 d6 vectors
 allowing for permulations of d1, d2, d3, d4
 */

#define CD7M_d7eucldistsq_byelem(v11,v12,v13,v14,v15,v16,v17,v21,v22,v23,v24,v25,v26,v27) \
fabs((v11-v21)*(v11-v21)+(v12-v22)*(v12-v22)+(v13-v23)*(v13-v23) + \
(v14-v24)*(v14-v24)+(v15-v25)*(v15-v25)+(v16-v26)*(v16-v26) + (v16-v26)*(v16-v26))

#define CD7M_d71234distsq(v1,v2) \
CD7M_min3(\
CD7M_min8(\
CD7M_d7eucldistsq_byelem(v1[0],v1[1],v1[2],v1[3],v1[4],v1[5],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[0],v1[1],v1[3],v1[2],v1[5],v1[4],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[0],v1[2],v1[1],v1[3],v1[4],v1[6],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[0],v1[2],v1[3],v1[1],v1[6],v1[4],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[0],v1[3],v1[1],v1[2],v1[5],v1[6],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[0],v1[3],v1[2],v1[1],v1[6],v1[5],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[2],v1[0],v1[1],v1[3],v1[6],v1[4],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[0],v1[1],v1[2],v1[6],v1[5],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6])),\
CD7M_min8(\
CD7M_d7eucldistsq_byelem(v1[1],v1[0],v1[2],v1[3],v1[5],v1[4],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[1],v1[0],v1[3],v1[2],v1[4],v1[5],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[1],v1[2],v1[0],v1[3],v1[5],v1[6],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[1],v1[2],v1[3],v1[0],v1[6],v1[5],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[1],v1[3],v1[0],v1[2],v1[4],v1[6],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[1],v1[3],v1[2],v1[0],v1[6],v1[4],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[2],v1[0],v1[3],v1[1],v1[4],v1[6],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[2],v1[1],v1[0],v1[3],v1[6],v1[5],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6])),\
CD7M_min8(\
CD7M_d7eucldistsq_byelem(v1[2],v1[1],v1[3],v1[0],v1[5],v1[6],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[2],v1[3],v1[0],v1[1],v1[4],v1[5],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[2],v1[3],v1[1],v1[0],v1[5],v1[4],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[0],v1[2],v1[1],v1[5],v1[6],v1[4],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[1],v1[0],v1[2],v1[6],v1[4],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[1],v1[2],v1[0],v1[4],v1[6],v1[5],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[2],v1[0],v1[1],v1[5],v1[4],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6]),\
CD7M_d7eucldistsq_byelem(v1[3],v1[2],v1[1],v1[0],v1[4],v1[5],v1[6],\
v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6])))

#define CD7M_d71234dist(v1,v2) sqrt(CD7M_d71234distsq(v1,v2))



/*     Compute the best distance between 2 G6 vectors
 allowing for permulations of g1, g2, g3 as
 well as sign changes
 */

#define CD7M_g456pmdistsq_byelem(v11,v12,v13,v14,v15,v16,v21,v22,v23,v24,v25,v26) \
fabs((v11-v21)*(v11-v21)+(v12-v22)*(v12-v22)+(v13-v23)*(v13-v23) + \
v14*v14+v24*v24 + v15*v15+v25*v25 + v16*v16+v26*v26 + \
2.*(CD7M_min3((fabs(v14*v24) - fabs(v15*v25) - fabs(v16*v26)), \
          ( - fabs(v14*v24) + fabs(v15*v25) - fabs(v16*v26)), \
          ( - fabs(v14*v24) - fabs(v15*v25) + fabs(v16*v26)))))


#define FASTER
/* #define FASTEST */
#ifdef FASTEST
#define CD7M_gtestdist(v1,v2) CD7M_gpmdist(v1,v2)
#define CD7M_gtestdistsq(v1,v2) CD7M_gpmdistsq(v1,v2)
#else
#ifdef FASTER
#define CD7M_gtestdist(v1,v2) CD7M_g456pmdist(v1,v2)
#define CD7M_gtestdistsq(v1,v2) CD7M_g456pmdistsq(v1,v2)
#else
#define CD7M_gtestdist(v1,v2) CD7M_g123pmdist(v1,v2)
#define CD7M_gtestdistsq(v1,v2) CD7M_g123pmdistsq(v1,v2)
#endif
#endif



static void d7cpyvn(int n, double src[], double dst[] ) {
    int i;
    for (i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}


static void imv7 (double v1[7], int m[49], double v2[7]) {
    int i, j;
    double sum;
    for (i = 0; i < 7; i++) {
        sum = 0.0;
        for(j=0; j < 7; j++) {
            sum = sum + ((double)m[7*i+j])*v1[j];
        }
        v2[i] = sum;
    }
}


static void rmv7 (double v1[7], double m[49], double v2[7]) {
    int i, j;
    double sum;
    for (i = 0; i < 7; i++) {
        sum = 0.0;
        for(j=0; j < 7; j++) {
            sum = sum + m[7*i+j]*v1[j];
        }
        v2[i] = sum;
    }
}




static void d7twoPminusI(double pg[7], double g[7], double gout[7]) {
    int i;
    for (i=0; i<7; i++){
        gout[i] = 2.*pg[i]-g[i];
    }
}


/*     Map a G6 vector onto the boundaries after
       applying the 24-way unfolding */

/* #define NCDIST_NO_OUTER_PASS */
#ifdef NCDIST_NO_OUTER_PASS
#define NREFL_OUTER_FULL 1
#define NREFL_OUTER_MIN 1
#else
#define NREFL_OUTER_FULL 24
#define NREFL_OUTER_MIN 6
#endif

static double d7bddist(double gvec[7],int bdnum) {
    
    double D7N_1x[7]={-19./sqrt(1106),-5./sqrt(1106),-12./sqrt(1106),-12./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_2x[7]={-12./sqrt(1106),-19./sqrt(1106),-5./sqrt(1106),-12./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_3x[7]={-12./sqrt(1106),-12./sqrt(1106),-19./sqrt(1106),-5./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_4x[7]={-3./sqrt(39),1./sqrt(39),1./sqrt(39),-3./sqrt(39),-1./sqrt(39),3./sqrt(39),3./sqrt(39)};
    double D7N_5x[7]={1./sqrt(39),-3./sqrt(39),-3./sqrt(39),1./sqrt(39),-1./sqrt(39),3./sqrt(39),3./sqrt(39)};
    double D7N_6x[7]={0,-1./2.,0,-1./2.,1./2.,0,1./2};
    double D7N_7x[7]={-3./sqrt(39),1./sqrt(39),-3./sqrt(39),1./sqrt(39),3./sqrt(39),-1./sqrt(39),3./sqrt(39)};
    double D7N_8x[7]={-2./4.,-2./4.,-1./4.,-1./4.,1./4.,1./4.,2./4.};
    double D7N_9x[7]={-3./sqrt(39),-3./sqrt(39),1./sqrt(39),1./sqrt(39),3./sqrt(39),3./sqrt(39),-1./sqrt(39)};
    
    
    if (bdnum < ND7BND) {
        
        switch(bdnum) {
            case(P_1): return d7dotprod(D7N_1x,gvec); break;
            case(P_2): return d7dotprod(D7N_2x,gvec); break;
            case(P_3): return d7dotprod(D7N_3x,gvec); break;
            case(P_4): return d7dotprod(D7N_4x,gvec); break;
            case(P_5): return d7dotprod(D7N_5x,gvec); break;
            case(P_6): return d7dotprod(D7N_6x,gvec); break;
            case(P_7): return d7dotprod(D7N_7x,gvec); break;
            case(P_8): return d7dotprod(D7N_8x,gvec); break;
            case(P_9): return d7dotprod(D7N_9x,gvec); break;
                
            default: return DBL_MAX; break;
        }
        
    }
    
    return DBL_MAX;
}




/* Get the minimum distance to the boundaries */

static double d7minbddist(double gvec[7]) {
    int ii;
    double dists[9];
    double minbd;
    dists[0] = fabs(d7bddist(gvec,P_1));
    dists[1] = fabs(d7bddist(gvec,P_2));
    dists[2] = fabs(d7bddist(gvec,P_3));
    dists[3] = fabs(d7bddist(gvec,P_4));
    dists[4] = fabs(d7bddist(gvec,P_5));
    dists[5] = fabs(d7bddist(gvec,P_6));
    dists[6] = fabs(d7bddist(gvec,P_7));
    dists[7] = fabs(d7bddist(gvec,P_8));
    dists[8] = fabs(d7bddist(gvec,P_9));
    
    minbd = dists[0];
    for (ii=1; ii<9; ii++) {
        minbd = CD7M_min(minbd,dists[ii]);
    }
    return minbd;
}



static void d7bdmaps(double gvec[7],
            double dists[ND7BND],
            int iord[ND7BND],
            double pgs[ND7BND][7],
            double rgs[ND7BND][7],
            double mpgs[ND7BND][7],
            double mvecs[ND7BND][7],
            double maxdist,
            int * ngood) {
    
    int jj, itemp, igap, idone;
    
    double D7N_1x[7]={-19./sqrt(1106),-5./sqrt(1106),-12./sqrt(1106),-12./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_2x[7]={-12./sqrt(1106),-19./sqrt(1106),-5./sqrt(1106),-12./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_3x[7]={-12./sqrt(1106),-12./sqrt(1106),-19./sqrt(1106),-5./sqrt(1106),12./sqrt(1106),12./sqrt(1106),12./sqrt(1106)};
    double D7N_4x[7]={-3./sqrt(39),1./sqrt(39),1./sqrt(39),-3./sqrt(39),-1./sqrt(39),3./sqrt(39),3./sqrt(39)};
    double D7N_5x[7]={1./sqrt(39),-3./sqrt(39),-3./sqrt(39),1./sqrt(39),-1./sqrt(39),3./sqrt(39),3./sqrt(39)};
    double D7N_6x[7]={0,-1./2.,0,-1./2.,1./2.,0,1./2};
    double D7N_7x[7]={-3./sqrt(39),1./sqrt(39),-3./sqrt(39),1./sqrt(39),3./sqrt(39),-1./sqrt(39),3./sqrt(39)};
    double D7N_8x[7]={-2./4.,-2./4.,-1./4.,-1./4.,1./4.,1./4.,2./4.};
    double D7N_9x[7]={-3./sqrt(39),-3./sqrt(39),1./sqrt(39),1./sqrt(39),3./sqrt(39),3./sqrt(39),-1./sqrt(39)};
    
    
    dists[P_1] = fabs(d7dotprod(D7N_1x,gvec));
    dists[P_2] = fabs(d7dotprod(D7N_2x,gvec));
    dists[P_3] = fabs(d7dotprod(D7N_3x,gvec));
    dists[P_4] = fabs(d7dotprod(D7N_4x,gvec));
    dists[P_5] = fabs(d7dotprod(D7N_5x,gvec));
    dists[P_6] = fabs(d7dotprod(D7N_6x,gvec));
    dists[P_7] = fabs(d7dotprod(D7N_7x,gvec));
    dists[P_8] = fabs(d7dotprod(D7N_8x,gvec));
    dists[P_9] = fabs(d7dotprod(D7N_9x,gvec));
                
    
    for (jj=0; jj < ND7BND; jj++) {
        iord[jj] = jj;
    }

    
    igap = ND7BND;
    while (igap > 1) {
        igap = igap/2;
        idone = 0;
        while (!idone) {
            idone = 1;
            for (jj=0; jj < ND7BND-igap; jj+=igap) {
                if (dists[iord[jj]] > dists[iord[jj+igap]]) {
                    idone = 0;
                    itemp = iord[jj];
                    iord[jj] = iord[jj+igap];
                    iord[jj+igap] = itemp;
                }
            }
        }
    }

    *ngood = ND7BND;
    for (jj = 0; jj < ND7BND; jj++ ) {
        rmv7(gvec, d7prj[jj], pgs[jj]);
        d7twoPminusI(pgs[jj],gvec,rgs[jj]);
        imv7(pgs[jj], D7MS[jj], mpgs[jj]);
        imv7(gvec, D7MS[jj],mvecs[jj]);
        if (dists[jj] > maxdist) (*ngood)--;
    }
    
}




#define DCUT 0.9995
#define fudge(d) DCUT*d


/*
 Compute tminimal distance between two Delaunay-reduced
 vectors in the Delaunay Cone following the embedding paths
 to the 9 boundaries
 */



double D7Dist_pass(double gvec1[7],double gvec2[7],double dist) {
    double dists1[ND7BND];
    double pgs1[ND7BND][7], rgs1[ND7BND][7], mpgs1[ND7BND][7], mvecs1[ND7BND][7];
    double dists2[ND7BND];
    double pgs2[ND7BND][7], rgs2[ND7BND][7], mpgs2[ND7BND][7], mvecs2[ND7BND][7];
    int iord1[ND7BND],iord2[ND7BND];
    double mindists1;
    double mindists2;
    int jx1, jx2;
    int j1,j2;
    int ngood1,ngood2;
    double maxdist;
    
    maxdist = fudge(dist);
    
    d7bdmaps(gvec1,dists1,iord1,pgs1,rgs1,mpgs1,mvecs1,maxdist,&ngood1);
    d7bdmaps(gvec2,dists2,iord2,pgs2,rgs2,mpgs2,mvecs2,maxdist,&ngood2);
    
    mindists1 = d7minbddist(gvec1);
    mindists2 = d7minbddist(gvec2);
    
    if (mindists1+mindists2 > dist) return dist;
    
    if (mindists1+mindists2 < maxdist) {
        for (jx1 = 0; jx1 < ngood1; jx1++) {
            double d1;
            j1 = iord1[jx1];
            d1 = dists1[j1];

            if (d1 < maxdist){
                dist = CD7M_min(dist,CD7M_d71234dist(gvec2,mpgs1[j1])+d1);
            }
        }
        for (jx2 = 0; jx2 < ngood2; jx2++) {
            double d2;
            j2 = iord2[jx2];
            d2 = dists2[j2];
            if (d2 < maxdist) {
                dist = CD7M_min(dist,(CD7M_d71234dist(gvec1,mpgs2[j2])+d2));
                
            }
        }
    }
    
    return dist;
}


/*
     Compute the minimal distance between two Niggli-reduced
     vectors in the Niggli Cone following the embedding paths
     to the 15 boundaries
 */





double D7Dist(double gvec1[7],double gvec2[7]) {
    double dist,dist1, dist2, distmin;
    dist1 = d7minbddist(gvec1);
    dist2 = d7minbddist(gvec2);
    distmin = CD7M_min(dist1,dist2);
    dist = CD7M_d71234dist(gvec1,gvec2);
    d7report_double("\n  Entered D7Dist d7dist = ",dist,", ");
    dist = D7Dist_pass(gvec1,gvec2,dist);
    d7report_double("dist1 = ",dist1,", ");
    d7report_double("dist2 = ",dist2,", ");
    return dist;
}



#endif /*D7DIST_H */


