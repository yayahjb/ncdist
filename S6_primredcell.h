
#ifdef __cplusplus

extern "C" {

#endif

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
    double * mat66);

#ifdef __cplusplus

}

#endif

