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
 **    filename:   stokes_q2p1_mf_operators_diag_def.c
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

void MatMultMF_Stokes_MixedFEM3d_diagB(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += FAC*(eta_gp*pow(dNudy[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudx[0],2));
  Y[1] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudy[0],2));
  Y[2] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudy[0],2) + 2.0*eta_gp*pow(dNudz[0],2));
  Y[3] += FAC*(eta_gp*pow(dNudy[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudx[1],2));
  Y[4] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudy[1],2));
  Y[5] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudy[1],2) + 2.0*eta_gp*pow(dNudz[1],2));
  Y[6] += FAC*(eta_gp*pow(dNudy[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudx[2],2));
  Y[7] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudy[2],2));
  Y[8] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudy[2],2) + 2.0*eta_gp*pow(dNudz[2],2));
  Y[9] += FAC*(eta_gp*pow(dNudy[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudx[3],2));
  Y[10] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudy[3],2));
  Y[11] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudy[3],2) + 2.0*eta_gp*pow(dNudz[3],2));
  Y[12] += FAC*(eta_gp*pow(dNudy[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudx[4],2));
  Y[13] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudy[4],2));
  Y[14] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudy[4],2) + 2.0*eta_gp*pow(dNudz[4],2));
  Y[15] += FAC*(eta_gp*pow(dNudy[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudx[5],2));
  Y[16] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudy[5],2));
  Y[17] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudy[5],2) + 2.0*eta_gp*pow(dNudz[5],2));
  Y[18] += FAC*(eta_gp*pow(dNudy[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudx[6],2));
  Y[19] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudy[6],2));
  Y[20] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudy[6],2) + 2.0*eta_gp*pow(dNudz[6],2));
  Y[21] += FAC*(eta_gp*pow(dNudy[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudx[7],2));
  Y[22] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudy[7],2));
  Y[23] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudy[7],2) + 2.0*eta_gp*pow(dNudz[7],2));
  Y[24] += FAC*(eta_gp*pow(dNudy[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudx[8],2));
  Y[25] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudy[8],2));
  Y[26] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudy[8],2) + 2.0*eta_gp*pow(dNudz[8],2));
  Y[27] += FAC*(eta_gp*pow(dNudy[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudx[9],2));
  Y[28] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudy[9],2));
  Y[29] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudy[9],2) + 2.0*eta_gp*pow(dNudz[9],2));
  Y[30] += FAC*(eta_gp*pow(dNudy[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudx[10],2));
  Y[31] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudy[10],2));
  Y[32] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudy[10],2) + 2.0*eta_gp*pow(dNudz[10],2));
  Y[33] += FAC*(eta_gp*pow(dNudy[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudx[11],2));
  Y[34] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudy[11],2));
  Y[35] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudy[11],2) + 2.0*eta_gp*pow(dNudz[11],2));
  Y[36] += FAC*(eta_gp*pow(dNudy[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudx[12],2));
  Y[37] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudy[12],2));
  Y[38] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudy[12],2) + 2.0*eta_gp*pow(dNudz[12],2));
  Y[39] += FAC*(eta_gp*pow(dNudy[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudx[13],2));
  Y[40] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudy[13],2));
  Y[41] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudy[13],2) + 2.0*eta_gp*pow(dNudz[13],2));
  Y[42] += FAC*(eta_gp*pow(dNudy[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudx[14],2));
  Y[43] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudy[14],2));
  Y[44] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudy[14],2) + 2.0*eta_gp*pow(dNudz[14],2));
  Y[45] += FAC*(eta_gp*pow(dNudy[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudx[15],2));
  Y[46] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudy[15],2));
  Y[47] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudy[15],2) + 2.0*eta_gp*pow(dNudz[15],2));
  Y[48] += FAC*(eta_gp*pow(dNudy[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudx[16],2));
  Y[49] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudy[16],2));
  Y[50] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudy[16],2) + 2.0*eta_gp*pow(dNudz[16],2));
  Y[51] += FAC*(eta_gp*pow(dNudy[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudx[17],2));
  Y[52] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudy[17],2));
  Y[53] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudy[17],2) + 2.0*eta_gp*pow(dNudz[17],2));
  Y[54] += FAC*(eta_gp*pow(dNudy[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudx[18],2));
  Y[55] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudy[18],2));
  Y[56] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudy[18],2) + 2.0*eta_gp*pow(dNudz[18],2));
  Y[57] += FAC*(eta_gp*pow(dNudy[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudx[19],2));
  Y[58] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudy[19],2));
  Y[59] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudy[19],2) + 2.0*eta_gp*pow(dNudz[19],2));
  Y[60] += FAC*(eta_gp*pow(dNudy[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudx[20],2));
  Y[61] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudy[20],2));
  Y[62] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudy[20],2) + 2.0*eta_gp*pow(dNudz[20],2));
  Y[63] += FAC*(eta_gp*pow(dNudy[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudx[21],2));
  Y[64] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudy[21],2));
  Y[65] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudy[21],2) + 2.0*eta_gp*pow(dNudz[21],2));
  Y[66] += FAC*(eta_gp*pow(dNudy[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudx[22],2));
  Y[67] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudy[22],2));
  Y[68] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudy[22],2) + 2.0*eta_gp*pow(dNudz[22],2));
  Y[69] += FAC*(eta_gp*pow(dNudy[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudx[23],2));
  Y[70] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudy[23],2));
  Y[71] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudy[23],2) + 2.0*eta_gp*pow(dNudz[23],2));
  Y[72] += FAC*(eta_gp*pow(dNudy[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudx[24],2));
  Y[73] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudy[24],2));
  Y[74] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudy[24],2) + 2.0*eta_gp*pow(dNudz[24],2));
  Y[75] += FAC*(eta_gp*pow(dNudy[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudx[25],2));
  Y[76] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudy[25],2));
  Y[77] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudy[25],2) + 2.0*eta_gp*pow(dNudz[25],2));
  Y[78] += FAC*(eta_gp*pow(dNudy[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudx[26],2));
  Y[79] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudy[26],2));
  Y[80] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudy[26],2) + 2.0*eta_gp*pow(dNudz[26],2));
  Y[81] += 0.0;
  Y[82] += 0.0;
  Y[83] += 0.0;
  Y[84] += 0.0;

// total operations = 1138
}
void MatMultMF_Stokes_MixedFEM3d_diagB11(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += FAC*(eta_gp*pow(dNudy[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudx[0],2));
  Y[1] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudy[0],2));
  Y[2] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudy[0],2) + 2.0*eta_gp*pow(dNudz[0],2));
  Y[3] += FAC*(eta_gp*pow(dNudy[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudx[1],2));
  Y[4] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudy[1],2));
  Y[5] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudy[1],2) + 2.0*eta_gp*pow(dNudz[1],2));
  Y[6] += FAC*(eta_gp*pow(dNudy[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudx[2],2));
  Y[7] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudy[2],2));
  Y[8] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudy[2],2) + 2.0*eta_gp*pow(dNudz[2],2));
  Y[9] += FAC*(eta_gp*pow(dNudy[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudx[3],2));
  Y[10] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudy[3],2));
  Y[11] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudy[3],2) + 2.0*eta_gp*pow(dNudz[3],2));
  Y[12] += FAC*(eta_gp*pow(dNudy[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudx[4],2));
  Y[13] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudy[4],2));
  Y[14] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudy[4],2) + 2.0*eta_gp*pow(dNudz[4],2));
  Y[15] += FAC*(eta_gp*pow(dNudy[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudx[5],2));
  Y[16] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudy[5],2));
  Y[17] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudy[5],2) + 2.0*eta_gp*pow(dNudz[5],2));
  Y[18] += FAC*(eta_gp*pow(dNudy[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudx[6],2));
  Y[19] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudy[6],2));
  Y[20] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudy[6],2) + 2.0*eta_gp*pow(dNudz[6],2));
  Y[21] += FAC*(eta_gp*pow(dNudy[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudx[7],2));
  Y[22] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudy[7],2));
  Y[23] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudy[7],2) + 2.0*eta_gp*pow(dNudz[7],2));
  Y[24] += FAC*(eta_gp*pow(dNudy[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudx[8],2));
  Y[25] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudy[8],2));
  Y[26] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudy[8],2) + 2.0*eta_gp*pow(dNudz[8],2));
  Y[27] += FAC*(eta_gp*pow(dNudy[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudx[9],2));
  Y[28] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudy[9],2));
  Y[29] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudy[9],2) + 2.0*eta_gp*pow(dNudz[9],2));
  Y[30] += FAC*(eta_gp*pow(dNudy[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudx[10],2));
  Y[31] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudy[10],2));
  Y[32] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudy[10],2) + 2.0*eta_gp*pow(dNudz[10],2));
  Y[33] += FAC*(eta_gp*pow(dNudy[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudx[11],2));
  Y[34] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudy[11],2));
  Y[35] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudy[11],2) + 2.0*eta_gp*pow(dNudz[11],2));
  Y[36] += FAC*(eta_gp*pow(dNudy[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudx[12],2));
  Y[37] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudy[12],2));
  Y[38] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudy[12],2) + 2.0*eta_gp*pow(dNudz[12],2));
  Y[39] += FAC*(eta_gp*pow(dNudy[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudx[13],2));
  Y[40] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudy[13],2));
  Y[41] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudy[13],2) + 2.0*eta_gp*pow(dNudz[13],2));
  Y[42] += FAC*(eta_gp*pow(dNudy[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudx[14],2));
  Y[43] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudy[14],2));
  Y[44] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudy[14],2) + 2.0*eta_gp*pow(dNudz[14],2));
  Y[45] += FAC*(eta_gp*pow(dNudy[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudx[15],2));
  Y[46] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudy[15],2));
  Y[47] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudy[15],2) + 2.0*eta_gp*pow(dNudz[15],2));
  Y[48] += FAC*(eta_gp*pow(dNudy[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudx[16],2));
  Y[49] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudy[16],2));
  Y[50] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudy[16],2) + 2.0*eta_gp*pow(dNudz[16],2));
  Y[51] += FAC*(eta_gp*pow(dNudy[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudx[17],2));
  Y[52] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudy[17],2));
  Y[53] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudy[17],2) + 2.0*eta_gp*pow(dNudz[17],2));
  Y[54] += FAC*(eta_gp*pow(dNudy[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudx[18],2));
  Y[55] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudy[18],2));
  Y[56] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudy[18],2) + 2.0*eta_gp*pow(dNudz[18],2));
  Y[57] += FAC*(eta_gp*pow(dNudy[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudx[19],2));
  Y[58] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudy[19],2));
  Y[59] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudy[19],2) + 2.0*eta_gp*pow(dNudz[19],2));
  Y[60] += FAC*(eta_gp*pow(dNudy[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudx[20],2));
  Y[61] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudy[20],2));
  Y[62] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudy[20],2) + 2.0*eta_gp*pow(dNudz[20],2));
  Y[63] += FAC*(eta_gp*pow(dNudy[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudx[21],2));
  Y[64] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudy[21],2));
  Y[65] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudy[21],2) + 2.0*eta_gp*pow(dNudz[21],2));
  Y[66] += FAC*(eta_gp*pow(dNudy[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudx[22],2));
  Y[67] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudy[22],2));
  Y[68] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudy[22],2) + 2.0*eta_gp*pow(dNudz[22],2));
  Y[69] += FAC*(eta_gp*pow(dNudy[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudx[23],2));
  Y[70] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudy[23],2));
  Y[71] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudy[23],2) + 2.0*eta_gp*pow(dNudz[23],2));
  Y[72] += FAC*(eta_gp*pow(dNudy[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudx[24],2));
  Y[73] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudy[24],2));
  Y[74] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudy[24],2) + 2.0*eta_gp*pow(dNudz[24],2));
  Y[75] += FAC*(eta_gp*pow(dNudy[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudx[25],2));
  Y[76] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudy[25],2));
  Y[77] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudy[25],2) + 2.0*eta_gp*pow(dNudz[25],2));
  Y[78] += FAC*(eta_gp*pow(dNudy[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudx[26],2));
  Y[79] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudy[26],2));
  Y[80] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudy[26],2) + 2.0*eta_gp*pow(dNudz[26],2));

// total operations = 1134
}
void MatMultMF_Stokes_MixedFEM3d_diagBuu(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += FAC*(eta_gp*pow(dNudy[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudx[0],2));
  Y[1] += FAC*(eta_gp*pow(dNudy[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudx[1],2));
  Y[2] += FAC*(eta_gp*pow(dNudy[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudx[2],2));
  Y[3] += FAC*(eta_gp*pow(dNudy[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudx[3],2));
  Y[4] += FAC*(eta_gp*pow(dNudy[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudx[4],2));
  Y[5] += FAC*(eta_gp*pow(dNudy[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudx[5],2));
  Y[6] += FAC*(eta_gp*pow(dNudy[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudx[6],2));
  Y[7] += FAC*(eta_gp*pow(dNudy[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudx[7],2));
  Y[8] += FAC*(eta_gp*pow(dNudy[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudx[8],2));
  Y[9] += FAC*(eta_gp*pow(dNudy[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudx[9],2));
  Y[10] += FAC*(eta_gp*pow(dNudy[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudx[10],2));
  Y[11] += FAC*(eta_gp*pow(dNudy[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudx[11],2));
  Y[12] += FAC*(eta_gp*pow(dNudy[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudx[12],2));
  Y[13] += FAC*(eta_gp*pow(dNudy[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudx[13],2));
  Y[14] += FAC*(eta_gp*pow(dNudy[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudx[14],2));
  Y[15] += FAC*(eta_gp*pow(dNudy[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudx[15],2));
  Y[16] += FAC*(eta_gp*pow(dNudy[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudx[16],2));
  Y[17] += FAC*(eta_gp*pow(dNudy[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudx[17],2));
  Y[18] += FAC*(eta_gp*pow(dNudy[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudx[18],2));
  Y[19] += FAC*(eta_gp*pow(dNudy[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudx[19],2));
  Y[20] += FAC*(eta_gp*pow(dNudy[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudx[20],2));
  Y[21] += FAC*(eta_gp*pow(dNudy[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudx[21],2));
  Y[22] += FAC*(eta_gp*pow(dNudy[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudx[22],2));
  Y[23] += FAC*(eta_gp*pow(dNudy[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudx[23],2));
  Y[24] += FAC*(eta_gp*pow(dNudy[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudx[24],2));
  Y[25] += FAC*(eta_gp*pow(dNudy[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudx[25],2));
  Y[26] += FAC*(eta_gp*pow(dNudy[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudx[26],2));

// total operations = 378
}
void MatMultMF_Stokes_MixedFEM3d_diagBvv(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudz[0],2) + 2.0*eta_gp*pow(dNudy[0],2));
  Y[1] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudz[1],2) + 2.0*eta_gp*pow(dNudy[1],2));
  Y[2] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudz[2],2) + 2.0*eta_gp*pow(dNudy[2],2));
  Y[3] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudz[3],2) + 2.0*eta_gp*pow(dNudy[3],2));
  Y[4] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudz[4],2) + 2.0*eta_gp*pow(dNudy[4],2));
  Y[5] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudz[5],2) + 2.0*eta_gp*pow(dNudy[5],2));
  Y[6] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudz[6],2) + 2.0*eta_gp*pow(dNudy[6],2));
  Y[7] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudz[7],2) + 2.0*eta_gp*pow(dNudy[7],2));
  Y[8] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudz[8],2) + 2.0*eta_gp*pow(dNudy[8],2));
  Y[9] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudz[9],2) + 2.0*eta_gp*pow(dNudy[9],2));
  Y[10] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudz[10],2) + 2.0*eta_gp*pow(dNudy[10],2));
  Y[11] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudz[11],2) + 2.0*eta_gp*pow(dNudy[11],2));
  Y[12] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudz[12],2) + 2.0*eta_gp*pow(dNudy[12],2));
  Y[13] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudz[13],2) + 2.0*eta_gp*pow(dNudy[13],2));
  Y[14] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudz[14],2) + 2.0*eta_gp*pow(dNudy[14],2));
  Y[15] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudz[15],2) + 2.0*eta_gp*pow(dNudy[15],2));
  Y[16] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudz[16],2) + 2.0*eta_gp*pow(dNudy[16],2));
  Y[17] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudz[17],2) + 2.0*eta_gp*pow(dNudy[17],2));
  Y[18] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudz[18],2) + 2.0*eta_gp*pow(dNudy[18],2));
  Y[19] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudz[19],2) + 2.0*eta_gp*pow(dNudy[19],2));
  Y[20] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudz[20],2) + 2.0*eta_gp*pow(dNudy[20],2));
  Y[21] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudz[21],2) + 2.0*eta_gp*pow(dNudy[21],2));
  Y[22] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudz[22],2) + 2.0*eta_gp*pow(dNudy[22],2));
  Y[23] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudz[23],2) + 2.0*eta_gp*pow(dNudy[23],2));
  Y[24] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudz[24],2) + 2.0*eta_gp*pow(dNudy[24],2));
  Y[25] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudz[25],2) + 2.0*eta_gp*pow(dNudy[25],2));
  Y[26] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudz[26],2) + 2.0*eta_gp*pow(dNudy[26],2));

// total operations = 378
}
void MatMultMF_Stokes_MixedFEM3d_diagBww(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += FAC*(eta_gp*pow(dNudx[0],2) + eta_gp*pow(dNudy[0],2) + 2.0*eta_gp*pow(dNudz[0],2));
  Y[1] += FAC*(eta_gp*pow(dNudx[1],2) + eta_gp*pow(dNudy[1],2) + 2.0*eta_gp*pow(dNudz[1],2));
  Y[2] += FAC*(eta_gp*pow(dNudx[2],2) + eta_gp*pow(dNudy[2],2) + 2.0*eta_gp*pow(dNudz[2],2));
  Y[3] += FAC*(eta_gp*pow(dNudx[3],2) + eta_gp*pow(dNudy[3],2) + 2.0*eta_gp*pow(dNudz[3],2));
  Y[4] += FAC*(eta_gp*pow(dNudx[4],2) + eta_gp*pow(dNudy[4],2) + 2.0*eta_gp*pow(dNudz[4],2));
  Y[5] += FAC*(eta_gp*pow(dNudx[5],2) + eta_gp*pow(dNudy[5],2) + 2.0*eta_gp*pow(dNudz[5],2));
  Y[6] += FAC*(eta_gp*pow(dNudx[6],2) + eta_gp*pow(dNudy[6],2) + 2.0*eta_gp*pow(dNudz[6],2));
  Y[7] += FAC*(eta_gp*pow(dNudx[7],2) + eta_gp*pow(dNudy[7],2) + 2.0*eta_gp*pow(dNudz[7],2));
  Y[8] += FAC*(eta_gp*pow(dNudx[8],2) + eta_gp*pow(dNudy[8],2) + 2.0*eta_gp*pow(dNudz[8],2));
  Y[9] += FAC*(eta_gp*pow(dNudx[9],2) + eta_gp*pow(dNudy[9],2) + 2.0*eta_gp*pow(dNudz[9],2));
  Y[10] += FAC*(eta_gp*pow(dNudx[10],2) + eta_gp*pow(dNudy[10],2) + 2.0*eta_gp*pow(dNudz[10],2));
  Y[11] += FAC*(eta_gp*pow(dNudx[11],2) + eta_gp*pow(dNudy[11],2) + 2.0*eta_gp*pow(dNudz[11],2));
  Y[12] += FAC*(eta_gp*pow(dNudx[12],2) + eta_gp*pow(dNudy[12],2) + 2.0*eta_gp*pow(dNudz[12],2));
  Y[13] += FAC*(eta_gp*pow(dNudx[13],2) + eta_gp*pow(dNudy[13],2) + 2.0*eta_gp*pow(dNudz[13],2));
  Y[14] += FAC*(eta_gp*pow(dNudx[14],2) + eta_gp*pow(dNudy[14],2) + 2.0*eta_gp*pow(dNudz[14],2));
  Y[15] += FAC*(eta_gp*pow(dNudx[15],2) + eta_gp*pow(dNudy[15],2) + 2.0*eta_gp*pow(dNudz[15],2));
  Y[16] += FAC*(eta_gp*pow(dNudx[16],2) + eta_gp*pow(dNudy[16],2) + 2.0*eta_gp*pow(dNudz[16],2));
  Y[17] += FAC*(eta_gp*pow(dNudx[17],2) + eta_gp*pow(dNudy[17],2) + 2.0*eta_gp*pow(dNudz[17],2));
  Y[18] += FAC*(eta_gp*pow(dNudx[18],2) + eta_gp*pow(dNudy[18],2) + 2.0*eta_gp*pow(dNudz[18],2));
  Y[19] += FAC*(eta_gp*pow(dNudx[19],2) + eta_gp*pow(dNudy[19],2) + 2.0*eta_gp*pow(dNudz[19],2));
  Y[20] += FAC*(eta_gp*pow(dNudx[20],2) + eta_gp*pow(dNudy[20],2) + 2.0*eta_gp*pow(dNudz[20],2));
  Y[21] += FAC*(eta_gp*pow(dNudx[21],2) + eta_gp*pow(dNudy[21],2) + 2.0*eta_gp*pow(dNudz[21],2));
  Y[22] += FAC*(eta_gp*pow(dNudx[22],2) + eta_gp*pow(dNudy[22],2) + 2.0*eta_gp*pow(dNudz[22],2));
  Y[23] += FAC*(eta_gp*pow(dNudx[23],2) + eta_gp*pow(dNudy[23],2) + 2.0*eta_gp*pow(dNudz[23],2));
  Y[24] += FAC*(eta_gp*pow(dNudx[24],2) + eta_gp*pow(dNudy[24],2) + 2.0*eta_gp*pow(dNudz[24],2));
  Y[25] += FAC*(eta_gp*pow(dNudx[25],2) + eta_gp*pow(dNudy[25],2) + 2.0*eta_gp*pow(dNudz[25],2));
  Y[26] += FAC*(eta_gp*pow(dNudx[26],2) + eta_gp*pow(dNudy[26],2) + 2.0*eta_gp*pow(dNudz[26],2));

// total operations = 378
}
void MatMultMF_Stokes_MixedFEM3d_diagA22(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  Y[0] += 0.0;
  Y[1] += 0.0;
  Y[2] += 0.0;
  Y[3] += 0.0;

// total operations = 4
}
