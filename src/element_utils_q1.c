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
 **    filename:   element_utils_q1.c
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

#include "element_utils_q1.h"

void P3D_ConstructNi_Q1_3D(PetscReal _xi[],PetscReal Ni[])
{
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];
  PetscReal zeta = _xi[2];

  Ni[0] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta ) * ( 1.0 - zeta );
  Ni[1] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta ) * ( 1.0 - zeta );
  Ni[2] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta ) * ( 1.0 - zeta );
  Ni[3] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta ) * ( 1.0 - zeta );

  Ni[4] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta ) * ( 1.0 + zeta );
  Ni[5] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta ) * ( 1.0 + zeta );
  Ni[6] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta ) * ( 1.0 + zeta );
  Ni[7] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta ) * ( 1.0 + zeta );
}

void P3D_ConstructGNi_Q1_3D(PetscReal _xi[],PetscReal GNi[3][8])
{
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];
  PetscReal zeta = _xi[2];

  GNi[0][0] = - 0.125 * ( 1.0 - eta ) * ( 1.0 - zeta );
  GNi[0][1] =   0.125 * ( 1.0 - eta ) * ( 1.0 - zeta );
  GNi[0][2] = - 0.125 * ( 1.0 + eta ) * ( 1.0 - zeta );
  GNi[0][3] =   0.125 * ( 1.0 + eta ) * ( 1.0 - zeta );

  GNi[0][4] = - 0.125 * ( 1.0 - eta ) * ( 1.0 + zeta );
  GNi[0][5] =   0.125 * ( 1.0 - eta ) * ( 1.0 + zeta );
  GNi[0][6] = - 0.125 * ( 1.0 + eta ) * ( 1.0 + zeta );
  GNi[0][7] =   0.125 * ( 1.0 + eta ) * ( 1.0 + zeta );
  //
  GNi[1][0] = - 0.125 * ( 1.0 - xi ) * ( 1.0 - zeta );
  GNi[1][1] = - 0.125 * ( 1.0 + xi ) * ( 1.0 - zeta );
  GNi[1][2] =   0.125 * ( 1.0 - xi ) * ( 1.0 - zeta );
  GNi[1][3] =   0.125 * ( 1.0 + xi ) * ( 1.0 - zeta );

  GNi[1][4] = - 0.125 * ( 1.0 - xi ) * ( 1.0 + zeta );
  GNi[1][5] = - 0.125 * ( 1.0 + xi ) * ( 1.0 + zeta );
  GNi[1][6] =   0.125 * ( 1.0 - xi ) * ( 1.0 + zeta );
  GNi[1][7] =   0.125 * ( 1.0 + xi ) * ( 1.0 + zeta );
  //
  GNi[2][0] = -0.125 * ( 1.0 - xi ) * ( 1.0 - eta );
  GNi[2][1] = -0.125 * ( 1.0 + xi ) * ( 1.0 - eta );
  GNi[2][2] = -0.125 * ( 1.0 - xi ) * ( 1.0 + eta );
  GNi[2][3] = -0.125 * ( 1.0 + xi ) * ( 1.0 + eta );

  GNi[2][4] = 0.125 * ( 1.0 - xi ) * ( 1.0 - eta );
  GNi[2][5] = 0.125 * ( 1.0 + xi ) * ( 1.0 - eta );
  GNi[2][6] = 0.125 * ( 1.0 - xi ) * ( 1.0 + eta );
  GNi[2][7] = 0.125 * ( 1.0 + xi ) * ( 1.0 + eta );
}


void P3D_evaluate_geometry_elementQ1(PetscInt nqp,PetscReal el_coords[8*3],PetscReal GNI[][3][8],
                                     PetscReal detJ[],
                                     PetscReal dNudx[][8],
                                     PetscReal dNudy[][8],
                                     PetscReal dNudz[][8] )
{
  PetscInt k,p;
  PetscReal t4, t6, t8, t10, t12, t14, t17;
  PetscReal J[3][3],iJ[3][3];

  for (p=0; p<nqp; p++) {
    //
    J[0][0] = J[0][1] = J[0][2] = 0.0;
    J[1][0] = J[1][1] = J[1][2] = 0.0;
    J[2][0] = J[2][1] = J[2][2] = 0.0;
    //
    //    memset(J[0],0,sizeof(PetscReal)*3);
    //    memset(J[1],0,sizeof(PetscReal)*3);
    //    memset(J[2],0,sizeof(PetscReal)*3);
    for (k=0; k<8; k++) {
      PetscReal xc = el_coords[3*k+0];
      PetscReal yc = el_coords[3*k+1];
      PetscReal zc = el_coords[3*k+2];

      J[0][0] += GNI[p][0][k] * xc ;
      J[0][1] += GNI[p][0][k] * yc ;
      J[0][2] += GNI[p][0][k] * zc ;

      J[1][0] += GNI[p][1][k] * xc ;
      J[1][1] += GNI[p][1][k] * yc ;
      J[1][2] += GNI[p][1][k] * zc ;

      J[2][0] += GNI[p][2][k] * xc ;
      J[2][1] += GNI[p][2][k] * yc ;
      J[2][2] += GNI[p][2][k] * zc ;
    }
    /* flops = [NQP*NPE] * 18 */

    detJ[p] = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
            - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
            + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
    /* flops = [NQP] * 14 */

    t4  = J[2][0] * J[0][1];
    t6  = J[2][0] * J[0][2];
    t8  = J[1][0] * J[0][1];
    t10 = J[1][0] * J[0][2];
    t12 = J[0][0] * J[1][1];
    t14 = J[0][0] * J[1][2]; // 6
    t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12

    iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
    iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
    iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
    iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
    iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
    iJ[1][2] = -(-t10 + t14) * t17;                            // 4
    iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
    iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
    iJ[2][2] = (-t8 + t12) * t17;                              // 3
    /* flops = [NQP] * 58 */

    /* shape function derivatives */
    for (k=0; k<8; k++) {
      dNudx[p][k] = iJ[0][0]*GNI[p][0][k] + iJ[0][1]*GNI[p][1][k] + iJ[0][2]*GNI[p][2][k];

      dNudy[p][k] = iJ[1][0]*GNI[p][0][k] + iJ[1][1]*GNI[p][1][k] + iJ[1][2]*GNI[p][2][k];

      dNudz[p][k] = iJ[2][0]*GNI[p][0][k] + iJ[2][1]*GNI[p][1][k] + iJ[2][2]*GNI[p][2][k];
    }
  }
  /* flops = [NQP*NPE] * 15 */
}

void P3D_evaluate_geometry_elementQ1_appliedQ2(PetscInt nqp,
                                     PetscReal detJ[],
                                     PetscReal GNIQ1[][3][8],
                                     PetscReal el_coords[27*3],
                                     PetscReal GNIQ2[][3][27],
                                     PetscReal dNudxQ2[][27],
                                     PetscReal dNudyQ2[][27],
                                     PetscReal dNudzQ2[][27] )
{
  PetscInt k,p;
  PetscReal t4, t6, t8, t10, t12, t14, t17;
  PetscReal J[3][3],iJ[3][3];
  PetscReal el_coordsQ1[8*3];

  for (k=0; k<3; k++) {
    el_coordsQ1[3*0+k] = el_coords[3*0+k];
    el_coordsQ1[3*1+k] = el_coords[3*2+k];
    el_coordsQ1[3*2+k] = el_coords[3*6+k];
    el_coordsQ1[3*3+k] = el_coords[3*8+k];

    el_coordsQ1[3*4+k] = el_coords[3*18+k];
    el_coordsQ1[3*5+k] = el_coords[3*20+k];
    el_coordsQ1[3*6+k] = el_coords[3*24+k];
    el_coordsQ1[3*7+k] = el_coords[3*26+k];
  }


  for (p=0; p<nqp; p++) {
    //
    J[0][0] = J[0][1] = J[0][2] = 0.0;
    J[1][0] = J[1][1] = J[1][2] = 0.0;
    J[2][0] = J[2][1] = J[2][2] = 0.0;
    //
    //    memset(J[0],0,sizeof(PetscReal)*3);
    //    memset(J[1],0,sizeof(PetscReal)*3);
    //    memset(J[2],0,sizeof(PetscReal)*3);
    for (k=0; k<8; k++) {
      PetscReal xc = el_coordsQ1[3*k+0];
      PetscReal yc = el_coordsQ1[3*k+1];
      PetscReal zc = el_coordsQ1[3*k+2];

      J[0][0] += GNIQ1[p][0][k] * xc ;
      J[0][1] += GNIQ1[p][0][k] * yc ;
      J[0][2] += GNIQ1[p][0][k] * zc ;

      J[1][0] += GNIQ1[p][1][k] * xc ;
      J[1][1] += GNIQ1[p][1][k] * yc ;
      J[1][2] += GNIQ1[p][1][k] * zc ;

      J[2][0] += GNIQ1[p][2][k] * xc ;
      J[2][1] += GNIQ1[p][2][k] * yc ;
      J[2][2] += GNIQ1[p][2][k] * zc ;
    }

    /* flops = [NQP*NPE] * 18 */

    detJ[p] = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
            - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
            + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
    /* flops = [NQP] * 14 */

    t4  = J[2][0] * J[0][1];
    t6  = J[2][0] * J[0][2];
    t8  = J[1][0] * J[0][1];
    t10 = J[1][0] * J[0][2];
    t12 = J[0][0] * J[1][1];
    t14 = J[0][0] * J[1][2]; // 6
    t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12

    iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
    iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
    iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
    iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
    iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
    iJ[1][2] = -(-t10 + t14) * t17;                            // 4
    iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
    iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
    iJ[2][2] = (-t8 + t12) * t17;                              // 3
    /* flops = [NQP] * 58 */

    /* shape function derivatives */
    for (k=0; k<27; k++) {
      dNudxQ2[p][k] = iJ[0][0]*GNIQ2[p][0][k] + iJ[0][1]*GNIQ2[p][1][k] + iJ[0][2]*GNIQ2[p][2][k];

      dNudyQ2[p][k] = iJ[1][0]*GNIQ2[p][0][k] + iJ[1][1]*GNIQ2[p][1][k] + iJ[1][2]*GNIQ2[p][2][k];

      dNudzQ2[p][k] = iJ[2][0]*GNIQ2[p][0][k] + iJ[2][1]*GNIQ2[p][1][k] + iJ[2][2]*GNIQ2[p][2][k];
    }
  }
  /* flops = [NQP*NPE] * 15 */

}

/* please refer to these papers */
/*

[1] On the nested refinement of quadrilateral and hexahedral finite elements and the affine approximation
  Shangyou Zhang
  Numer. Math. (2004) 98: 559–579

[2] Numerical integration with Taylor truncations for the quadrilateral and hexahedral finite elements
  Shangyou Zhang
  Journal of Computational and Applied Mathematics 205 (2007) 325 – 342

[3] A STABLE AFFINE-APPROXIMATE FINITE ELEMENT METHOD
  K. ARUNAKIRINATHAR AND B. D. REDDY
  SIAM J. NUMER. ANAL.
  Vol. 40, No. 1, pp. 180–197

*/
void P3D_evaluate_geometry_affine_appliedQ2(PetscInt nqp,
                                     PetscReal detJ[],
                                     PetscReal GNIQ1[][3][8],
                                     PetscReal el_coords[27*3],
                                     PetscReal GNIQ2[][3][27],
                                     PetscReal dNudxQ2[][27],
                                     PetscReal dNudyQ2[][27],
                                     PetscReal dNudzQ2[][27] )
{
  PetscInt k,p;
  PetscReal t4, t6, t8, t10, t12, t14, t17;
  PetscReal J[3][3],iJ[3][3],detJ_p;
  PetscReal a100[3],a010[3],a001[3];
  PetscReal v1[3],v2[3],v3[3],v4[3],v5[3],v6[3],v7[3],v8[3];


  for (k=0; k<3; k++) {
    v1[k] = el_coords[3*0+k];
    v2[k] = el_coords[3*2+k];
    v3[k] = el_coords[3*6+k];
    v4[k] = el_coords[3*8+k];

    v5[k] = el_coords[3*18+k];
    v6[k] = el_coords[3*20+k];
    v7[k] = el_coords[3*24+k];
    v8[k] = el_coords[3*26+k];
  }

/*

 MAP = a000 + a100x + a010y + a001z

 a000 =  0.125*v1 + 0.125*v2 + 0.125*v3 + 0.125*v4 + 0.125*v5 + 0.125*v6 + 0.125*v7 + 0.125*v8
 a100 =  0.125*v3 + 0.125*v4 + 0.125*v7 + 0.125*v8 - 0.125*v1 - 0.125*v2 - 0.125*v5 - 0.125*v6
 a010 =  0.125*v2 + 0.125*v3 + 0.125*v6 + 0.125*v7 - 0.125*v1 - 0.125*v4 - 0.125*v5 - 0.125*v8
 a001 =  0.125*v5 + 0.125*v6 + 0.125*v7 + 0.125*v8 - 0.125*v1 - 0.125*v2 - 0.125*v3 - 0.125*v4

*/

  for (k=0; k<3; k++) {
    a100[k] =  0.125*(v2[k] + v4[k] + v6[k] + v8[k] - v1[k] - v3[k] - v5[k] - v7[k]);
    a010[k] =  0.125*(v3[k] + v4[k] + v7[k] + v8[k] - v1[k] - v2[k] - v5[k] - v6[k]);
    a001[k] =  0.125*(v5[k] + v6[k] + v7[k] + v8[k] - v1[k] - v2[k] - v3[k] - v4[k]);
  }

  J[0][0] = 1.0*a100[0] ;
  J[0][1] = 1.0*a010[0] ;
  J[0][2] = 1.0*a001[0] ;

  J[1][0] = 1.0*a100[1] ;
  J[1][1] = 1.0*a010[1] ;
  J[1][2] = 1.0*a001[1] ;

  J[2][0] = 1.0*a100[2] ;
  J[2][1] = 1.0*a010[2] ;
  J[2][2] = 1.0*a001[2] ;

  detJ_p = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) // a
         - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
         + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]); // c
  /* flops =  14 */
  //printf("detJ = %1.4e \n",detJ_p);

  for (p=0; p<nqp; p++) {
    detJ[p] = detJ_p;
  }


  t4  = J[2][0] * J[0][1];
  t6  = J[2][0] * J[0][2];
  t8  = J[1][0] * J[0][1];
  t10 = J[1][0] * J[0][2];
  t12 = J[0][0] * J[1][1];
  t14 = J[0][0] * J[1][2]; // 6
  t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1]);  // 12

  iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17;  // 4
  iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17; // 5
  iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17;  // 4
  iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17;// 6
  iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17;                // 4
  iJ[1][2] = -(-t10 + t14) * t17;                            // 4
  iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17; // 5
  iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17;               // 5
  iJ[2][2] = (-t8 + t12) * t17;                              // 3
  /* flops = 58 */


  for (p=0; p<nqp; p++) {
    /* shape function derivatives */
    for (k=0; k<27; k++) {
      dNudxQ2[p][k] = iJ[0][0]*GNIQ2[p][0][k] + iJ[0][1]*GNIQ2[p][1][k] + iJ[0][2]*GNIQ2[p][2][k];

      dNudyQ2[p][k] = iJ[1][0]*GNIQ2[p][0][k] + iJ[1][1]*GNIQ2[p][1][k] + iJ[1][2]*GNIQ2[p][2][k];

      dNudzQ2[p][k] = iJ[2][0]*GNIQ2[p][0][k] + iJ[2][1]*GNIQ2[p][1][k] + iJ[2][2]*GNIQ2[p][2][k];
    }
  }
  /* flops = [NQP*NPE] * 15 */

}


void P3D_evaluate_geometry_affine2_appliedQ2(PetscInt nqp,
                                        PetscReal detJ[],
                                        PetscReal GNIQ1[][3][8],
                                        PetscReal el_coords[27*3],
                                        PetscReal GNIQ2[][3][27],
                                        PetscReal dNudxQ2[][27],
                                        PetscReal dNudyQ2[][27],
                                        PetscReal dNudzQ2[][27] )
{
  PetscInt k,p;
  PetscReal J[3],iJ[3],detJ_p;
  PetscReal a100[3],a010[3],a001[3];
  PetscReal v1[3],v2[3],v3[3],v4[3],v5[3],v6[3],v7[3],v8[3];

  for (k=0; k<3; k++) {
    v1[k] = el_coords[3*0+k];
    v2[k] = el_coords[3*2+k];
    v3[k] = el_coords[3*6+k];
    v4[k] = el_coords[3*8+k];

    v5[k] = el_coords[3*18+k];
    v6[k] = el_coords[3*20+k];
    v7[k] = el_coords[3*24+k];
    v8[k] = el_coords[3*26+k];
  }

  /*

   MAP = a000 + a100x + a010y + a001z

   a000 =  0.125*v1 + 0.125*v2 + 0.125*v3 + 0.125*v4 + 0.125*v5 + 0.125*v6 + 0.125*v7 + 0.125*v8
   a100 =  0.125*v2 + 0.125*v4 + 0.125*v6 + 0.125*v8 - 0.125*v1 - 0.125*v3 - 0.125*v5 - 0.125*v7
   a010 =  0.125*v3 + 0.125*v4 + 0.125*v7 + 0.125*v8 - 0.125*v1 - 0.125*v2 - 0.125*v5 - 0.125*v6
   a001 =  0.125*v5 + 0.125*v6 + 0.125*v7 + 0.125*v8 - 0.125*v1 - 0.125*v2 - 0.125*v3 - 0.125*v4

   */

  for (k=0; k<3; k++) {
    a100[k] =  0.125*(v2[k] + v4[k] + v6[k] + v8[k] - v1[k] - v3[k] - v5[k] - v7[k]);
    a010[k] =  0.125*(v3[k] + v4[k] + v7[k] + v8[k] - v1[k] - v2[k] - v5[k] - v6[k]);
    a001[k] =  0.125*(v5[k] + v6[k] + v7[k] + v8[k] - v1[k] - v2[k] - v3[k] - v4[k]);
  }


  J[0] = 1.0*a100[0] ;
  J[1] = 1.0*a010[1] ;
  J[2] = 1.0*a001[2] ;

  detJ_p = J[0]*J[1]*J[2];
  /* flops =  2 */

  for (p=0; p<nqp; p++) {
    detJ[p] = detJ_p;
  }

  iJ[0] = 1.0/J[0];
  iJ[1] = 1.0/J[1];
  iJ[2] = 1.0/J[2];


  for (p=0; p<nqp; p++) {
    /* shape function derivatives */
    for (k=0; k<27; k++) {
      dNudxQ2[p][k] = iJ[0]*GNIQ2[p][0][k];

      dNudyQ2[p][k] = iJ[1]*GNIQ2[p][1][k];

      dNudzQ2[p][k] = iJ[2]*GNIQ2[p][2][k];
    }
  }

}


void P3D_ConstructNi_Q1_2D(PetscReal _xi[],PetscReal Ni[])
{
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];

  Ni[0] = 0.25 * ( 1.0 - xi ) * ( 1.0 - eta );
  Ni[1] = 0.25 * ( 1.0 + xi ) * ( 1.0 - eta );
  Ni[2] = 0.25 * ( 1.0 - xi ) * ( 1.0 + eta );
  Ni[3] = 0.25 * ( 1.0 + xi ) * ( 1.0 + eta );
}

void P3D_ConstructGNi_Q1_2D(PetscReal _xi[],PetscReal GNix[],PetscReal GNiy[])
{
  PetscReal xi   = _xi[0];
  PetscReal eta  = _xi[1];

  GNix[0] = - 0.25 * ( 1.0 - eta );
  GNix[1] =   0.25 * ( 1.0 - eta );
  GNix[2] = - 0.25 * ( 1.0 + eta );
  GNix[3] =   0.25 * ( 1.0 + eta );
  //
  GNiy[0] = - 0.25 * ( 1.0 - xi );
  GNiy[1] = - 0.25 * ( 1.0 + xi );
  GNiy[2] =   0.25 * ( 1.0 - xi );
  GNiy[3] =   0.25 * ( 1.0 + xi );
}

void P3D_evaluate_geometry_elementQ1_2D(PetscReal el_coords[],PetscReal GNIx[],PetscReal GNIy[],
                                     PetscReal *detJ,PetscReal dNudx[],PetscReal dNudy[])
{
  PetscInt k;
  PetscReal dJ,J[2][2],iJ[2][2];

    J[0][0] = J[0][1] = 0.0;
    J[1][0] = J[1][1] = 0.0;

    for (k=0; k<4; k++) {
        PetscReal xc = el_coords[2*k+0];
        PetscReal yc = el_coords[2*k+1];

        J[0][0] += GNIx[k] * xc ;
        J[0][1] += GNIx[k] * yc ;

        J[1][0] += GNIy[k] * xc ;
        J[1][1] += GNIy[k] * yc ;
    }
    dJ = (J[0][0]*J[1][1] - J[0][1]*J[1][0]);

    iJ[0][0] =  J[1][1]/dJ;
    iJ[0][1] = -J[0][1]/dJ;
    iJ[1][0] = -J[1][0]/dJ;
    iJ[1][1] =  J[0][0]/dJ;
    *detJ = dJ;

    /* shape function derivatives */
    for (k=0; k<4; k++) {
        dNudx[k] = iJ[0][0]*GNIx[k] + iJ[0][1]*GNIy[k];
        dNudy[k] = iJ[1][0]*GNIx[k] + iJ[1][1]*GNIy[k];
    }
}


