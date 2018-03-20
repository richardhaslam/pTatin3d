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
 **    filename:   stokes_q2p1_mf_operators_def_rolled.c
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

#if defined(PTAT3D_ELEMENT_MF_OPTIMIZED)
#include <petsclog.h>           /* for PetscLogFlops() */

static inline void MatMultMF_Stokes_MixedFEM3d_B(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu,ip,k;
  double    p_gp = 0.0;
  double    strain_rate_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    tau_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    div_gp = 0.0;
  double    te,tmp;


  // pressure at gauss point //
  for (ip=0; ip<4; ip++) {
    p_gp += Np[ip]*P[ip];
  } // ops = 4.[2]

  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[0] += Ux[iu]*dNudx[iu];
    strain_rate_gp[1] += Uy[iu]*dNudy[iu];
    strain_rate_gp[2] += Uz[iu]*dNudz[iu];
    strain_rate_gp[3] += Ux[iu]*dNudy[iu] + Uy[iu]*dNudx[iu];
    strain_rate_gp[4] += Ux[iu]*dNudz[iu] + Uz[iu]*dNudx[iu];
    strain_rate_gp[5] += Uy[iu]*dNudz[iu] + Uz[iu]*dNudy[iu];
  } // ops = 27.[18]

  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  te = 2.0 * eta_gp;
  tau_gp[0] = te*strain_rate_gp[0];
  tau_gp[1] = te*strain_rate_gp[1];
  tau_gp[2] = te*strain_rate_gp[2];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // ops = 1 + 6 = 7

  // divergence at gauss point //
  for (k=0; k<3; k++) {
    div_gp += strain_rate_gp[k];
  } // ops = 3

  // y = A11.u + A12.p at gauss point //
  for (iu=0; iu<27; iu++) {
    int idx = 3*iu;
    Y[idx]   += FAC*(dNudx[iu]*tau_gp[0] + dNudy[iu]*tau_gp[3] + dNudz[iu]*tau_gp[4] - dNudx[iu]*p_gp);
    Y[++idx] += FAC*(dNudx[iu]*tau_gp[3] + dNudy[iu]*tau_gp[1] + dNudz[iu]*tau_gp[5] - dNudy[iu]*p_gp);
    Y[++idx] += FAC*(dNudx[iu]*tau_gp[4] + dNudy[iu]*tau_gp[5] + dNudz[iu]*tau_gp[2] - dNudz[iu]*p_gp);
  } // ops = 27.[27]

  // y = A21.u at gauss point //
  tmp = -FAC * div_gp;
  for (ip=0; ip<4; ip++) {
    int idx = 81+ip;
    Y[idx] += tmp*Np[ip];
  } // ops = 2 + 4.[2]
  // total operations = 8 + 490 + 7 + 729 + 10 = 1244
  PetscLogFlops(1244);
}
static inline void MatMultMF_Stokes_MixedFEM3d_B11(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu;
  double    strain_rate_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    tau_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    te;


  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[0] += Ux[iu]*dNudx[iu];
    strain_rate_gp[1] += Uy[iu]*dNudy[iu];
    strain_rate_gp[2] += Uz[iu]*dNudz[iu];
    strain_rate_gp[3] += Ux[iu]*dNudy[iu] + Uy[iu]*dNudx[iu];
    strain_rate_gp[4] += Ux[iu]*dNudz[iu] + Uz[iu]*dNudx[iu];
    strain_rate_gp[5] += Uy[iu]*dNudz[iu] + Uz[iu]*dNudy[iu];
  } // ops = 27.[18]

  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  te = 2.0 * eta_gp;
  tau_gp[0] = te*strain_rate_gp[0];
  tau_gp[1] = te*strain_rate_gp[1];
  tau_gp[2] = te*strain_rate_gp[2];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // ops = 1 + 6 = 7

  // y = A11.u at gauss point //
  for (iu=0; iu<27; iu++) {
    int idx = 3*iu;
    Y[idx]   += (dNudx[iu]*tau_gp[0] + dNudy[iu]*tau_gp[3] + dNudz[iu]*tau_gp[4]);
    Y[++idx] += (dNudx[iu]*tau_gp[3] + dNudy[iu]*tau_gp[1] + dNudz[iu]*tau_gp[5]);
    Y[++idx] += (dNudx[iu]*tau_gp[4] + dNudy[iu]*tau_gp[5] + dNudz[iu]*tau_gp[2]);
  } // ops = 27.[18]

// total operations = 486 + 7 + 486 = 979
  PetscLogFlops(979);
}
static inline void MatMultMF_Stokes_MixedFEM3d_Buu(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu;
  double    strain_rate_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    tau_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};


  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[0] += Ux[iu]*dNudx[iu];
    strain_rate_gp[3] += Ux[iu]*dNudy[iu];
    strain_rate_gp[4] += Ux[iu]*dNudz[iu];
  }

  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 2.0*eta_gp*strain_rate_gp[0];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];

  // y = A11.{u,0,0} at gauss point //
  for (iu=0; iu<27; iu++) {
    Y[iu]  += FAC*(dNudx[iu]*tau_gp[0] + dNudy[iu]*tau_gp[3] + dNudz[iu]*tau_gp[4]);
  }
}
static inline void MatMultMF_Stokes_MixedFEM3d_Bvv(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu;
  double    strain_rate_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    tau_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};


  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[1] += Uy[iu]*dNudy[iu];
    strain_rate_gp[3] += Uy[iu]*dNudx[iu];
    strain_rate_gp[5] += Uy[iu]*dNudz[iu];
  }

  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[1] = 2.0*eta_gp*strain_rate_gp[1];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[5] = eta_gp*strain_rate_gp[5];

  // y = A11.{0,v,0} at gauss point //
  for (iu=0; iu<27; iu++) {
    Y[iu]  += FAC*(dNudx[iu]*tau_gp[3] + dNudy[iu]*tau_gp[1] + dNudz[iu]*tau_gp[5]);
  }
}
static inline void MatMultMF_Stokes_MixedFEM3d_Bww(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu;
  double    strain_rate_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};
  double    tau_gp[] = {0.0,0.0,0.0,0.0,0.0,0.0};


  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[2] += Uz[iu]*dNudz[iu];
    strain_rate_gp[4] += Uz[iu]*dNudx[iu];
    strain_rate_gp[5] += Uz[iu]*dNudy[iu];
  }

  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[2] = 2.0*eta_gp*strain_rate_gp[2];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];

  // y = A11.{0,0,w} at gauss point //
  for (iu=0; iu<27; iu++) {
    Y[iu]  += FAC*(dNudx[iu]*tau_gp[4] + dNudy[iu]*tau_gp[5] + dNudz[iu]*tau_gp[2]);
  }
}
static inline void MatMultMF_Stokes_MixedFEM3d_A12(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu,ip;
  double    p_gp = 0.0,tmp;


  // pressure at gauss point //
  for (ip=0; ip<4; ip++) {
    p_gp += Np[ip]*P[ip];
  }

  // y = A12.p at gauss point //
  tmp = -FAC * p_gp;
  for (iu=0; iu<27; iu++) {
    int idx = 3*iu;
    Y[  idx] += tmp*dNudx[iu];
    Y[++idx] += tmp*dNudy[iu];
    Y[++idx] += tmp*dNudz[iu];
  }
}
static inline void MatMultMF_Stokes_MixedFEM3d_A21(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd   = 3;
  const int ntens = 6;
  int       iu,ip;
  double    strain_rate_gp[] = {0.0,0.0,0.0};
  double    div_gp,tmp;


  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  for (iu=0; iu<27; iu++) {
    strain_rate_gp[0] += Ux[iu]*dNudx[iu];
    strain_rate_gp[1] += Uy[iu]*dNudy[iu];
    strain_rate_gp[2] += Uz[iu]*dNudz[iu];
  }

  // divergence at gauss point //
  div_gp = strain_rate_gp[0] + strain_rate_gp[1] + strain_rate_gp[2];

  // y = A21.u at gauss point //
  tmp = -FAC*div_gp;
  for (ip=0; ip<4; ip++) {
    Y[ip] += tmp*Np[ip];
  }
}
static inline void MatMultMF_Stokes_MixedFEM3d_A22(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  // y = A22.p at gauss point //
  // do nothing! //
  /*
  Y[0] = 0;
  Y[1] = 0;
  Y[2] = 0;
  Y[3] = 0;
  */
}
#endif
