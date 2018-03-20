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
 **    filename:   compute_analytic_growth_rate_3d.c
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
/*
 Plot using gnuplot

 sets contours from 2-40 in increments of 2

 gnuplot> set contour base
 gnuplot> set cntrparam levels incremental 2,2,40
 gnuplot> splot "l.dat" u 1:2:4 w l

*/


#include "stdio.h"
#include "math.h"

int main(void)
{
  double exx,eyy,ezz;
  double q,q3d,R,kx,kz,k,eta_l,eta_m,nu;
  double alpha,beta;
  int i,j;

  eta_l = 100.0;
  eta_m = 1.0;

  // R = 1
  exx = -1.0;
  ezz = -0.0;

  // R = 0.75
  exx = -1.0;
  ezz = -0.333333333333333;

  // R = 1.5
  //exx = -1.0;
  //ezz = 0.333333333333333;


  eyy = -(exx + ezz);

  nu = eta_l/eta_m;
  R = - exx / eyy;

  printf("# R = %1.4e \n",R);
  printf("# nu = %1.4e \n",nu);

  for (j=1; j<100; j++) {
    for (i=1; i<100; i++) {
      double lambda_x = i * (0.6);
      double lambda_z = j * (0.6);

      kx = (2.0 * M_PI) / lambda_x;
      kz = (2.0 * M_PI) / lambda_z;

      k = sqrt(kx*kx + kz*kz);

      alpha = -4.0 * (1.0 - 1.0/nu) * k;
      beta  = 2.0 * k * (1.0 - 1.0/(nu*nu));
      beta -= (1.0 + 1.0/nu)*(1.0 + 1.0/nu) * exp(k);
      beta += (1.0 - 1.0/nu)*(1.0 - 1.0/nu) * exp(-k);
      q = alpha / beta;

      q3d = eyy * 0.5 * q * ( kz*kz*(R-1.0)/(k*k) - kx*kx*R/(k*k) - 1.0 );

      printf("%1.4e %1.4e %1.4e %1.4e\n", lambda_x,lambda_z,q3d,-q3d/eyy);
    }
    printf("\n");
  }

  return 0;
}

