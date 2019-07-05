#if _MSC_VER
#pragma once
#endif // _MSC_VER

#ifndef RANDTOOLS_H
#define RANDTOOLS_H

#include "LRL_Cell.h"
#include "rhrand.h"
#include "G6.h"


void initrn( int& iseed, int& indx, int& jndx, int& kndx, double buffer[56] );

void SetSeed( const int n );

double randfg( int& iseed );

G6 GenRandG6( void );

G6 genRandG6V() ;

LRL_Cell GenRandG6ReducedLRL_Cell( void );
LRL_Cell GenRandG6DeloneReducedLRL_Cell( void );

G6 GenRandG6LRL_Cell( void );

bool GenRandReducedArray(double* a, double aRan[6]);
bool GenRandDeloneReducedArray(double* a, double aRan[6]);

#endif

