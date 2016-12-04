                          ncdist -- Niggli Cone Distance utilities
                        Herbert J. Bernstein and Lawrence C. Andrews
                                Version 1.0.0 03 Dec 2016

                                    yayahjb@gmail.com


ncdist is a package of two utilities built on the Niggli-Cone cell distance discussed in 
[Andrews, L. C.,  Bernstein, H. J. (2014). "The geometry of Niggli reduction: BGAOL -- 
embedding Niggli reduction and analysis of boundaries". Journal of applied crystallography, 
47(1), 346-359.] which provides a sensitive and reliable metric for the distance between
crystallographic unit cells.  A shell command line program ncdist is provided as well as
a shared library rcpp_ncdist for use of ncdist from R.

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

THE COMMAND LINE PROGRAM NCDIST

The program ncdist accepts two sets of 7 command line arguments:

  ncdist lat1 a1 b1 c1 alpha1 beta1 gamma1 lat2 a2 b2 c2 alpha2 beta2 gamma2

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

A call to

ncdist P 100 100 100 90 90 90 P 99 99 99 89 89 89

will produce

Primitive Reduced Probe Cell: 100 100 100 90 90 90
Volume :1e+06
Reciprocal of Primitive Probe Cell: 0.01 0.01 0.01 90 90 90
Volume of Reciprocal Cell: 1e-06
V7 linearized and scaled: 92.582 92.582 92.582 92.582 92.582 92.582  92.582
raw G6 vector: 10000 10000 10000 1.22465e-12 1.22465e-12 1.22465e-12

Primitive Reduced Probe Cell: 99 99 99 89 89 89
Volume :969861
Reciprocal of Primitive Probe Cell: 0.010104 0.010104 0.010104 90.9828 90.9828 90.9828
Volume of Reciprocal Cell: 1.03108e-06
V7 linearized and scaled: 91.6562 91.6562 91.6562 91.6287 91.6287 91.6287  91.6424
raw G6 vector: 9801 9801 9801 342.102 342.102 342.102

dprim1: [10000, 10000, 10000, 1.22465e-12, 1.22465e-12, 10000]
dprim2: [9801, 9801, 9801, 342.102, 342.102, 9801]
2.6182

Which says the NCdist distance between these two cells if 2.6182 Angstroms, computed
as 0.1 times the square root of the G6 distance in Angstroms^2.  At present the
extra information above the distance is provided to facilitate debugging and will
be made optional in a future version.

THE RCPP INTERFACE TO NCDIST

The library rcpp_ncdist is provided as an interface between R and ncdist.  Depending on
the type of system used, the name of the file may be a .dylib, a .dll, a .so and may or may
not be prefixed by lib.  Under linux the name is librcpp_ncdist.so.  Under OS X the name
is librcpp_ncdist.dylib.

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

