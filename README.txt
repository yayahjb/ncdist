                          ncdist -- Niggli Cone Distance utilities
                                 With S6 and D7 Support
                        Herbert J. Bernstein and Lawrence C. Andrews
                                Version 1.1.0 07 Jun 2019

                                    yayahjb@gmail.com


ncdist is a package of two utilities built on the Niggli-Cone cell distance discussed in 
[Andrews, L. C.,  Bernstein, H. J. (2014). "The geometry of Niggli reduction: BGAOL -- 
embedding Niggli reduction and analysis of boundaries". J. Appl. Cryst. 
47(1), 346 -- 359.] which provides a sensitive and reliable metric for the distance between
crystallographic unit cells.  Subsequent work has shown that similar distances can be
computed more efficiently based on Selling reduction which is the basis for Delaunay
reduction [Andrews, L. C., Bernstein, H. J. and Sauter, N. K. (2019). "Selling reduction 
versus Niggli reduction for crystallographic lattices". Acta Cryst. A75(1), 115 -- 120.]
[Andrews, L. C., Bernstein, H. J. and Sauter, N. K. (2019). "A space for lattice 
representation and clustering". Acta Cryst. A75(3), 593 -- 599.  Four shell command line 
programs, ncdist, d7dist, s6dist_app, and cs6dist_app are provided as well as
a shared libraries, rcpp_ncdist (librcpp_ncdist.so), rcpp_d7dist (librcpp_d7dist.so),
rcpp_s6dist (librcpp_s6dist.so),  and rcpp_cs6dist (librcpp_cs6dist.so),  for use from R.

The 1.1.0 release of 7 Jun 2019 has been refactored to use code from L. C. Andrews
Lattice Representation Library (https://github.com/duck10/LatticeRepLib.git) with
changes to adapt to the conventions of the ncdist code.

Work on the upgrade to version 1.1.0 has been supported in part by funding from the US 
Department of Energy Offices of Biological and Environmental Research and of Basic Energy 
Sciences (grant No. DE-AC02-98CH10886 ; grant No. E-SC0012704); U.S. National Institutes 
of Health (grant No. P41RR012408; grant No. P41GM103473; grant No. P41GM111244; and grant 
No. R01GM117126 ); Dectris, Ltd.

LICENSE

#####################################################################
#                                                                    #
# YOU MAY REDISTRIBUTE THE ncdist PACKAGE UNDER THE TERMS OF THE GPL #
#                                                                    #
# ALTERNATIVELY YOU MAY REDISTRIBUTE THE ncdist API UNDER THE TERMS  #
# OF THE LGPL                                                        #
#                                                                    #
######################################################################

########################### GPL NOTICES ##############################
#                                                                    #
# This program is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU General Public License as     #
# published by the Free Software Foundation; either version 2 of     #
# (the License, or (at your option) any later version.               #
#                                                                    #
# This program is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      #
# GNU General Public License for more details.                       #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with this program; if not, write to the Free Software        #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA           #
# 02111-1307  USA                                                    #
#                                                                    #
######################################################################

######################### LGPL NOTICES ###############################
#                                                                    #
# This library is free software; you can redistribute it and/or      #
# modify it under the terms of the GNU Lesser General Public         #
# License as published by the Free Software Foundation; either       #
# version 2.1 of the License, or (at your option) any later version. #
#                                                                    #
# This library is distributed in the hope that it will be useful,    #
# but WITHOUT ANY WARRANTY; without even the implied warranty of     #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# Lesser General Public License for more details.                    #
#                                                                    #
# You should have received a copy of the GNU Lesser General Public   #
# License along with this library; if not, write to the Free         #
# Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    #
# MA  02110-1301  USA                                                #
#                                                                    #
######################################################################

THE COMMAND LINE PROGRAMS NCDIST

The programs ncdist, d7dist, s6dist_app and cs6dist_app accept two sets of 7 command 
line arguments:

  ncdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2
  d7dist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2
  s6dist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2
  cs6dist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2

lat1 and lat2 are centering symbols
  P (Primitive)
  A (a-centered)
  B (b-centered)
  C (c-centered)
  F (all-faces-centered)
  I (body-centered)
  R (rhombohedral as hexagonal)
  H (hexagonal)
a1, b1, c1, a2, b2, c2 are cell edge length in Angstroms
alpha1, beta1, gamma1, alpha2, beta2, gamma2 are cell angles in degrees

A timed call to ncdist on a 12-core AMD threadripper will produce

time ncdist P 100 100 100 90 90 90 P 99 99 99 89 89 89

Primitive Reduced Probe Cell: 100 100 100 90 90 90
Delaunay Primitive Reduced Probe Cell: 100 100 100 90 90 90
Volume :1e+06
Reciprocal of Primitive Probe Cell: 0.01 0.01 0.01 90 90 90
Volume of Reciprocal Cell: 1e-06
V7 linearized and scaled: 92.582 92.582 92.582 92.582 92.582 92.582  92.582
raw G6 vector: 10000 10000 10000 1.45473e-11 1.45473e-11 1.45473e-11

Primitive Reduced Probe Cell: 99 99 99 89 89 89
Delaunay Primitive Reduced Probe Cell: 99 99 138.78 134.5 90 91
Volume :969861
Reciprocal of Primitive Probe Cell: 0.010104 0.010104 0.010104 90.9828 90.9828 90.9828
Volume of Reciprocal Cell: 1.03108e-06
V7 linearized and scaled: 91.6562 91.6562 91.6562 91.6287 91.6287 91.6287  91.6424
raw G6 vector: 9801 9801 9801 342.102 342.102 342.102

dprim1: [10000, 10000, 10000, 0, 0, 0]
dprim2: [9801, 9801, 9801, 342.102, 342.102, 4.68822e-310]
2.6182

real	0m0.006s
user	0m0.028s
sys	0m0.000s

Similar calls to d7dist, s6dist_app, and cs6dist_app produce

time ./d7dist P 100 100 100 90 90 90 P 99 99 99 89 89 89

dprim1: [10000, 10000, 10000, 30000, 20000, 20000, 20000]
dprim2: [9801, 9801, 19259.9, 19259.9, 9801, 29060.9, 20000]
4.0324

real	0m0.009s
user	0m0.044s
sys	0m0.004s

time ./s6dist_app P 100 100 100 90 90 90 P 99 99 99 89 89 89

dprim1: [0, 0, 0, -10000, -10000, -10000]
dprim2: [-9629.95, -9629.95, 0, -9629.95, -171.051, 0]
rawdist: 663.379
2.57561

real	0m0.005s
user	0m0.000s
sys	0m0.004s

time ./cs6dist_app P 100 100 100 90 90 90 P 99 99 99 89 89 89

Primitive Reduced Probe LRL_Cell: 100 100 100 90 90 90
Delaunay Primitive Reduced Probe LRL_Cell: 100 100 100 90 90 90
Selling Primitive Reduced Probe S6 LRL_Cell: 0 0 0 -10000 -10000 -10000
Volume :1e+06
Reciprocal of Primitive Probe LRL_Cell: 0.0001 0.0001 0.0001 -8.96147e-05 -8.96147e-05 -8.96147e-05
Volume of Reciprocal LRL_Cell: 1e-06
V7 linearized and scaled: 92.582 92.582 92.582 92.582 92.582 92.582  92.582
raw G6 vector for S6: 10000 10000 10000 1.22465e-12 1.22465e-12 1.22465e-12

Primitive Reduced Probe LRL_Cell: 99 99 99 89 89 89
Delaunay Primitive Reduced Probe LRL_Cell: 99 99 138.78 134.5 90 91
Selling Primitive Reduced Probe S6 LRL_Cell: -9629.95 -9629.95 0 -9629.95 -171.051 0
Volume :969861
Reciprocal of Primitive Probe LRL_Cell: 0.000102092 0.000102092 0.000102092 -0.000202632 -0.000202632 -0.000202632
Volume of Reciprocal LRL_Cell: 1.03108e-06
V7 linearized and scaled: 91.6562 91.6562 91.6562 91.6287 91.6287 91.6287  91.6424
raw G6 vector for S6: 19259.9 9801 19259.9 -19259.9 -19259.9 1.68257e-12

dprim1: [0, 0, 0, -10000, -10000, -10000]
dprim2: [-9629.95, -9629.95, 0, -9629.95, -171.051, 0]
rawdist: 663.379
2.57561

real	0m0.018s
user	0m0.280s
sys	0m0.000s


Which says the NCdist distance between these two cells if 2.6182 Angstroms, computed
as 0.1 times the square root of the G6 distance in Angstroms^2.  At present the
extra information above the distance is provided to facilitate debugging and will
be made optional in a future version.

THE RCPP INTERFACE TO NCDIST

The libraries rcpp_ncdist, rcpp_d7dist, rcpp_s6dist, and rcpp_cs6dist are provided as 
an interface between R and ncdist.  Depending on the type of system used, the name of 
the file may be a .dylib, a .dll, a .so and may or may not be prefixed by lib.  Under 
linux the names are librcpp_ncdist.so, librcpp_d7dist.so, librcpp_s6dist.so, and 
librcpp_cs6dist.so.  Under OS X the names are usually librcpp_ncdist.dylib, 
librcpp_ncdist.dylib,librcpp_s6dist.dylib, and librcpp_cs6dist.dylib

Assuming the library is named /usr/local/lib/librcpp_ncdist.so , in order to compute 
a cell distance in R, the commands are:

> library("Rcpp")
> dyn.load("/usr/local/lib/librcpp_ncdist.so")
> .Call("rcpp_ncdist", "P", 100, 100, 100, 90, 90, 90, "P", 99, 99, 99, 89, 89, 89)

which will return

[1] 2.618197

INSTALLATION

Preliminaries:  You need a development system with cmake, C, C++, Fortran, R, Rcpp, 
RcppParallel and RcppArmadillo installed.

Download the package from https://github.com/yayahjb/ncdist.git

If you downloaded as a zip, unpack the kit.  In any case if you have the source in
ncdist

cd ncdist/build
cmake -DCMAKE_INSTALL_PREFIX=<prefix> ..
make all
make install

will install the kit in <prefix>, specifically, it will 

install ncdist in <prefix>/bin/ncdist
install librcpp_ncdist in <prefix>/lib/librcpp_ncdist.<shared_library_extension>
install license.txt in <prefix>/share/ncdist/license.txt

In most cases you will want <prefix> to be $CCP4

