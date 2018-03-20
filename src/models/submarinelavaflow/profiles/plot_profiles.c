/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2012
 **        Dave A. May [dave.may@erdw.ethz.ch]
 **        Institute of Geophysics
 **        ETH Zürich
 **        Sonneggstrasse 5
 **        CH-8092 Zürich
 **        Switzerland
 **
 **    project:    pTatin3d
 **    filename:   plot_profiles.c
 **
 **
 **    pTatin3d is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, either version 3 of the License,
 **    or (at your option) any later version.
 **
 **    pTatin3d is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with pTatin3d. If not, see <http://www.gnu.org/licenses/>.
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main(int nargs,char *args[])
{
  int i,n;
  double x,x0,x1,dx;
  double phi_0,phi_f,Ts,Te,phi_max;
  double eta0,gamma;
  double T,phi,sy,eta;

  x0 = 0.0;
  x1 = 1.0;

  n = 10000;
  dx = (x1-x0)/((double)(n-1));

  for (i=0; i<n; i++) {
    x = x0 + i*dx;

    if (x <= 0.5) {
      T = 1100.0;
    } else {
      T = 0.0;
    }
    T = 1100.0 - i * 1100.0/((double)(n-1));

    phi_0 = 0.0;
    phi_max = 0.68;
    Ts = 600.0;
    Te = 1100.0;
    phi_f = phi_max - phi_0;
    phi = phi_0 + phi_f * (Te - T)/(Te - Ts);

    if (phi >= phi_max) { phi = phi_max - 1.0e-10; }

    eta0 = 50.0;
    gamma = -0.04;
    eta = eta0 * pow( 1.0 - phi/phi_max, -2.5 ) * exp( -gamma * (Te - T) );
    //if (eta > 1.0e5) { eta = 1.0e5; }

    sy = pow( 10.0, 11.59 - 0.0089*T );

    printf("%1.4e %1.4e %1.4e %1.4e %1.4e\n",x,T,phi,eta,sy);

  }


  return 0;
}

