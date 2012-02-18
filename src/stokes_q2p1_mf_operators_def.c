

inline void MatMultMF_Stokes_MixedFEM3d_B(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = Np[0]*P[0] + Np[1]*P[1] + Np[2]*P[2] + Np[3]*P[3];
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = Ux[0]*dNudx[0] + Ux[10]*dNudx[10] + Ux[11]*dNudx[11] + Ux[12]*dNudx[12] + Ux[13]*dNudx[13] + Ux[14]*dNudx[14] + Ux[15]*dNudx[15] + Ux[16]*dNudx[16] + Ux[17]*dNudx[17] + Ux[18]*dNudx[18] + Ux[19]*dNudx[19] + Ux[1]*dNudx[1] + Ux[20]*dNudx[20] + Ux[21]*dNudx[21] + Ux[22]*dNudx[22] + Ux[23]*dNudx[23] + Ux[24]*dNudx[24] + Ux[25]*dNudx[25] + Ux[26]*dNudx[26] + Ux[2]*dNudx[2] + Ux[3]*dNudx[3] + Ux[4]*dNudx[4] + Ux[5]*dNudx[5] + Ux[6]*dNudx[6] + Ux[7]*dNudx[7] + Ux[8]*dNudx[8] + Ux[9]*dNudx[9];
  strain_rate_gp[1] = Uy[0]*dNudy[0] + Uy[10]*dNudy[10] + Uy[11]*dNudy[11] + Uy[12]*dNudy[12] + Uy[13]*dNudy[13] + Uy[14]*dNudy[14] + Uy[15]*dNudy[15] + Uy[16]*dNudy[16] + Uy[17]*dNudy[17] + Uy[18]*dNudy[18] + Uy[19]*dNudy[19] + Uy[1]*dNudy[1] + Uy[20]*dNudy[20] + Uy[21]*dNudy[21] + Uy[22]*dNudy[22] + Uy[23]*dNudy[23] + Uy[24]*dNudy[24] + Uy[25]*dNudy[25] + Uy[26]*dNudy[26] + Uy[2]*dNudy[2] + Uy[3]*dNudy[3] + Uy[4]*dNudy[4] + Uy[5]*dNudy[5] + Uy[6]*dNudy[6] + Uy[7]*dNudy[7] + Uy[8]*dNudy[8] + Uy[9]*dNudy[9];
  strain_rate_gp[2] = Uz[0]*dNudz[0] + Uz[10]*dNudz[10] + Uz[11]*dNudz[11] + Uz[12]*dNudz[12] + Uz[13]*dNudz[13] + Uz[14]*dNudz[14] + Uz[15]*dNudz[15] + Uz[16]*dNudz[16] + Uz[17]*dNudz[17] + Uz[18]*dNudz[18] + Uz[19]*dNudz[19] + Uz[1]*dNudz[1] + Uz[20]*dNudz[20] + Uz[21]*dNudz[21] + Uz[22]*dNudz[22] + Uz[23]*dNudz[23] + Uz[24]*dNudz[24] + Uz[25]*dNudz[25] + Uz[26]*dNudz[26] + Uz[2]*dNudz[2] + Uz[3]*dNudz[3] + Uz[4]*dNudz[4] + Uz[5]*dNudz[5] + Uz[6]*dNudz[6] + Uz[7]*dNudz[7] + Uz[8]*dNudz[8] + Uz[9]*dNudz[9];
  strain_rate_gp[3] = Ux[0]*dNudy[0] + Ux[10]*dNudy[10] + Ux[11]*dNudy[11] + Ux[12]*dNudy[12] + Ux[13]*dNudy[13] + Ux[14]*dNudy[14] + Ux[15]*dNudy[15] + Ux[16]*dNudy[16] + Ux[17]*dNudy[17] + Ux[18]*dNudy[18] + Ux[19]*dNudy[19] + Ux[1]*dNudy[1] + Ux[20]*dNudy[20] + Ux[21]*dNudy[21] + Ux[22]*dNudy[22] + Ux[23]*dNudy[23] + Ux[24]*dNudy[24] + Ux[25]*dNudy[25] + Ux[26]*dNudy[26] + Ux[2]*dNudy[2] + Ux[3]*dNudy[3] + Ux[4]*dNudy[4] + Ux[5]*dNudy[5] + Ux[6]*dNudy[6] + Ux[7]*dNudy[7] + Ux[8]*dNudy[8] + Ux[9]*dNudy[9] + Uy[0]*dNudx[0] + Uy[10]*dNudx[10] + Uy[11]*dNudx[11] + Uy[12]*dNudx[12] + Uy[13]*dNudx[13] + Uy[14]*dNudx[14] + Uy[15]*dNudx[15] + Uy[16]*dNudx[16] + Uy[17]*dNudx[17] + Uy[18]*dNudx[18] + Uy[19]*dNudx[19] + Uy[1]*dNudx[1] + Uy[20]*dNudx[20] + Uy[21]*dNudx[21] + Uy[22]*dNudx[22] + Uy[23]*dNudx[23] + Uy[24]*dNudx[24] + Uy[25]*dNudx[25] + Uy[26]*dNudx[26] + Uy[2]*dNudx[2] + Uy[3]*dNudx[3] + Uy[4]*dNudx[4] + Uy[5]*dNudx[5] + Uy[6]*dNudx[6] + Uy[7]*dNudx[7] + Uy[8]*dNudx[8] + Uy[9]*dNudx[9];
  strain_rate_gp[4] = Ux[0]*dNudz[0] + Ux[10]*dNudz[10] + Ux[11]*dNudz[11] + Ux[12]*dNudz[12] + Ux[13]*dNudz[13] + Ux[14]*dNudz[14] + Ux[15]*dNudz[15] + Ux[16]*dNudz[16] + Ux[17]*dNudz[17] + Ux[18]*dNudz[18] + Ux[19]*dNudz[19] + Ux[1]*dNudz[1] + Ux[20]*dNudz[20] + Ux[21]*dNudz[21] + Ux[22]*dNudz[22] + Ux[23]*dNudz[23] + Ux[24]*dNudz[24] + Ux[25]*dNudz[25] + Ux[26]*dNudz[26] + Ux[2]*dNudz[2] + Ux[3]*dNudz[3] + Ux[4]*dNudz[4] + Ux[5]*dNudz[5] + Ux[6]*dNudz[6] + Ux[7]*dNudz[7] + Ux[8]*dNudz[8] + Ux[9]*dNudz[9] + Uz[0]*dNudx[0] + Uz[10]*dNudx[10] + Uz[11]*dNudx[11] + Uz[12]*dNudx[12] + Uz[13]*dNudx[13] + Uz[14]*dNudx[14] + Uz[15]*dNudx[15] + Uz[16]*dNudx[16] + Uz[17]*dNudx[17] + Uz[18]*dNudx[18] + Uz[19]*dNudx[19] + Uz[1]*dNudx[1] + Uz[20]*dNudx[20] + Uz[21]*dNudx[21] + Uz[22]*dNudx[22] + Uz[23]*dNudx[23] + Uz[24]*dNudx[24] + Uz[25]*dNudx[25] + Uz[26]*dNudx[26] + Uz[2]*dNudx[2] + Uz[3]*dNudx[3] + Uz[4]*dNudx[4] + Uz[5]*dNudx[5] + Uz[6]*dNudx[6] + Uz[7]*dNudx[7] + Uz[8]*dNudx[8] + Uz[9]*dNudx[9];
  strain_rate_gp[5] = Uy[0]*dNudz[0] + Uy[10]*dNudz[10] + Uy[11]*dNudz[11] + Uy[12]*dNudz[12] + Uy[13]*dNudz[13] + Uy[14]*dNudz[14] + Uy[15]*dNudz[15] + Uy[16]*dNudz[16] + Uy[17]*dNudz[17] + Uy[18]*dNudz[18] + Uy[19]*dNudz[19] + Uy[1]*dNudz[1] + Uy[20]*dNudz[20] + Uy[21]*dNudz[21] + Uy[22]*dNudz[22] + Uy[23]*dNudz[23] + Uy[24]*dNudz[24] + Uy[25]*dNudz[25] + Uy[26]*dNudz[26] + Uy[2]*dNudz[2] + Uy[3]*dNudz[3] + Uy[4]*dNudz[4] + Uy[5]*dNudz[5] + Uy[6]*dNudz[6] + Uy[7]*dNudz[7] + Uy[8]*dNudz[8] + Uy[9]*dNudz[9] + Uz[0]*dNudy[0] + Uz[10]*dNudy[10] + Uz[11]*dNudy[11] + Uz[12]*dNudy[12] + Uz[13]*dNudy[13] + Uz[14]*dNudy[14] + Uz[15]*dNudy[15] + Uz[16]*dNudy[16] + Uz[17]*dNudy[17] + Uz[18]*dNudy[18] + Uz[19]*dNudy[19] + Uz[1]*dNudy[1] + Uz[20]*dNudy[20] + Uz[21]*dNudy[21] + Uz[22]*dNudy[22] + Uz[23]*dNudy[23] + Uz[24]*dNudy[24] + Uz[25]*dNudy[25] + Uz[26]*dNudy[26] + Uz[2]*dNudy[2] + Uz[3]*dNudy[3] + Uz[4]*dNudy[4] + Uz[5]*dNudy[5] + Uz[6]*dNudy[6] + Uz[7]*dNudy[7] + Uz[8]*dNudy[8] + Uz[9]*dNudy[9];
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 2.0*eta_gp*strain_rate_gp[0];
  tau_gp[1] = 2.0*eta_gp*strain_rate_gp[1];
  tau_gp[2] = 2.0*eta_gp*strain_rate_gp[2];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[0] + strain_rate_gp[1] + strain_rate_gp[2];
  // y = A21.u at gauss point //
  Y[0] += dNudx[0]*tau_gp[0] + dNudy[0]*tau_gp[3] + dNudz[0]*tau_gp[4] - dNudx[0]*p_gp;
  Y[1] += dNudx[0]*tau_gp[3] + dNudy[0]*tau_gp[1] + dNudz[0]*tau_gp[5] - dNudy[0]*p_gp;
  Y[2] += dNudx[0]*tau_gp[4] + dNudy[0]*tau_gp[5] + dNudz[0]*tau_gp[2] - dNudz[0]*p_gp;
  Y[3] += dNudx[1]*tau_gp[0] + dNudy[1]*tau_gp[3] + dNudz[1]*tau_gp[4] - dNudx[1]*p_gp;
  Y[4] += dNudx[1]*tau_gp[3] + dNudy[1]*tau_gp[1] + dNudz[1]*tau_gp[5] - dNudy[1]*p_gp;
  Y[5] += dNudx[1]*tau_gp[4] + dNudy[1]*tau_gp[5] + dNudz[1]*tau_gp[2] - dNudz[1]*p_gp;
  Y[6] += dNudx[2]*tau_gp[0] + dNudy[2]*tau_gp[3] + dNudz[2]*tau_gp[4] - dNudx[2]*p_gp;
  Y[7] += dNudx[2]*tau_gp[3] + dNudy[2]*tau_gp[1] + dNudz[2]*tau_gp[5] - dNudy[2]*p_gp;
  Y[8] += dNudx[2]*tau_gp[4] + dNudy[2]*tau_gp[5] + dNudz[2]*tau_gp[2] - dNudz[2]*p_gp;
  Y[9] += dNudx[3]*tau_gp[0] + dNudy[3]*tau_gp[3] + dNudz[3]*tau_gp[4] - dNudx[3]*p_gp;
  Y[10] += dNudx[3]*tau_gp[3] + dNudy[3]*tau_gp[1] + dNudz[3]*tau_gp[5] - dNudy[3]*p_gp;
  Y[11] += dNudx[3]*tau_gp[4] + dNudy[3]*tau_gp[5] + dNudz[3]*tau_gp[2] - dNudz[3]*p_gp;
  Y[12] += dNudx[4]*tau_gp[0] + dNudy[4]*tau_gp[3] + dNudz[4]*tau_gp[4] - dNudx[4]*p_gp;
  Y[13] += dNudx[4]*tau_gp[3] + dNudy[4]*tau_gp[1] + dNudz[4]*tau_gp[5] - dNudy[4]*p_gp;
  Y[14] += dNudx[4]*tau_gp[4] + dNudy[4]*tau_gp[5] + dNudz[4]*tau_gp[2] - dNudz[4]*p_gp;
  Y[15] += dNudx[5]*tau_gp[0] + dNudy[5]*tau_gp[3] + dNudz[5]*tau_gp[4] - dNudx[5]*p_gp;
  Y[16] += dNudx[5]*tau_gp[3] + dNudy[5]*tau_gp[1] + dNudz[5]*tau_gp[5] - dNudy[5]*p_gp;
  Y[17] += dNudx[5]*tau_gp[4] + dNudy[5]*tau_gp[5] + dNudz[5]*tau_gp[2] - dNudz[5]*p_gp;
  Y[18] += dNudx[6]*tau_gp[0] + dNudy[6]*tau_gp[3] + dNudz[6]*tau_gp[4] - dNudx[6]*p_gp;
  Y[19] += dNudx[6]*tau_gp[3] + dNudy[6]*tau_gp[1] + dNudz[6]*tau_gp[5] - dNudy[6]*p_gp;
  Y[20] += dNudx[6]*tau_gp[4] + dNudy[6]*tau_gp[5] + dNudz[6]*tau_gp[2] - dNudz[6]*p_gp;
  Y[21] += dNudx[7]*tau_gp[0] + dNudy[7]*tau_gp[3] + dNudz[7]*tau_gp[4] - dNudx[7]*p_gp;
  Y[22] += dNudx[7]*tau_gp[3] + dNudy[7]*tau_gp[1] + dNudz[7]*tau_gp[5] - dNudy[7]*p_gp;
  Y[23] += dNudx[7]*tau_gp[4] + dNudy[7]*tau_gp[5] + dNudz[7]*tau_gp[2] - dNudz[7]*p_gp;
  Y[24] += dNudx[8]*tau_gp[0] + dNudy[8]*tau_gp[3] + dNudz[8]*tau_gp[4] - dNudx[8]*p_gp;
  Y[25] += dNudx[8]*tau_gp[3] + dNudy[8]*tau_gp[1] + dNudz[8]*tau_gp[5] - dNudy[8]*p_gp;
  Y[26] += dNudx[8]*tau_gp[4] + dNudy[8]*tau_gp[5] + dNudz[8]*tau_gp[2] - dNudz[8]*p_gp;
  Y[27] += dNudx[9]*tau_gp[0] + dNudy[9]*tau_gp[3] + dNudz[9]*tau_gp[4] - dNudx[9]*p_gp;
  Y[28] += dNudx[9]*tau_gp[3] + dNudy[9]*tau_gp[1] + dNudz[9]*tau_gp[5] - dNudy[9]*p_gp;
  Y[29] += dNudx[9]*tau_gp[4] + dNudy[9]*tau_gp[5] + dNudz[9]*tau_gp[2] - dNudz[9]*p_gp;
  Y[30] += dNudx[10]*tau_gp[0] + dNudy[10]*tau_gp[3] + dNudz[10]*tau_gp[4] - dNudx[10]*p_gp;
  Y[31] += dNudx[10]*tau_gp[3] + dNudy[10]*tau_gp[1] + dNudz[10]*tau_gp[5] - dNudy[10]*p_gp;
  Y[32] += dNudx[10]*tau_gp[4] + dNudy[10]*tau_gp[5] + dNudz[10]*tau_gp[2] - dNudz[10]*p_gp;
  Y[33] += dNudx[11]*tau_gp[0] + dNudy[11]*tau_gp[3] + dNudz[11]*tau_gp[4] - dNudx[11]*p_gp;
  Y[34] += dNudx[11]*tau_gp[3] + dNudy[11]*tau_gp[1] + dNudz[11]*tau_gp[5] - dNudy[11]*p_gp;
  Y[35] += dNudx[11]*tau_gp[4] + dNudy[11]*tau_gp[5] + dNudz[11]*tau_gp[2] - dNudz[11]*p_gp;
  Y[36] += dNudx[12]*tau_gp[0] + dNudy[12]*tau_gp[3] + dNudz[12]*tau_gp[4] - dNudx[12]*p_gp;
  Y[37] += dNudx[12]*tau_gp[3] + dNudy[12]*tau_gp[1] + dNudz[12]*tau_gp[5] - dNudy[12]*p_gp;
  Y[38] += dNudx[12]*tau_gp[4] + dNudy[12]*tau_gp[5] + dNudz[12]*tau_gp[2] - dNudz[12]*p_gp;
  Y[39] += dNudx[13]*tau_gp[0] + dNudy[13]*tau_gp[3] + dNudz[13]*tau_gp[4] - dNudx[13]*p_gp;
  Y[40] += dNudx[13]*tau_gp[3] + dNudy[13]*tau_gp[1] + dNudz[13]*tau_gp[5] - dNudy[13]*p_gp;
  Y[41] += dNudx[13]*tau_gp[4] + dNudy[13]*tau_gp[5] + dNudz[13]*tau_gp[2] - dNudz[13]*p_gp;
  Y[42] += dNudx[14]*tau_gp[0] + dNudy[14]*tau_gp[3] + dNudz[14]*tau_gp[4] - dNudx[14]*p_gp;
  Y[43] += dNudx[14]*tau_gp[3] + dNudy[14]*tau_gp[1] + dNudz[14]*tau_gp[5] - dNudy[14]*p_gp;
  Y[44] += dNudx[14]*tau_gp[4] + dNudy[14]*tau_gp[5] + dNudz[14]*tau_gp[2] - dNudz[14]*p_gp;
  Y[45] += dNudx[15]*tau_gp[0] + dNudy[15]*tau_gp[3] + dNudz[15]*tau_gp[4] - dNudx[15]*p_gp;
  Y[46] += dNudx[15]*tau_gp[3] + dNudy[15]*tau_gp[1] + dNudz[15]*tau_gp[5] - dNudy[15]*p_gp;
  Y[47] += dNudx[15]*tau_gp[4] + dNudy[15]*tau_gp[5] + dNudz[15]*tau_gp[2] - dNudz[15]*p_gp;
  Y[48] += dNudx[16]*tau_gp[0] + dNudy[16]*tau_gp[3] + dNudz[16]*tau_gp[4] - dNudx[16]*p_gp;
  Y[49] += dNudx[16]*tau_gp[3] + dNudy[16]*tau_gp[1] + dNudz[16]*tau_gp[5] - dNudy[16]*p_gp;
  Y[50] += dNudx[16]*tau_gp[4] + dNudy[16]*tau_gp[5] + dNudz[16]*tau_gp[2] - dNudz[16]*p_gp;
  Y[51] += dNudx[17]*tau_gp[0] + dNudy[17]*tau_gp[3] + dNudz[17]*tau_gp[4] - dNudx[17]*p_gp;
  Y[52] += dNudx[17]*tau_gp[3] + dNudy[17]*tau_gp[1] + dNudz[17]*tau_gp[5] - dNudy[17]*p_gp;
  Y[53] += dNudx[17]*tau_gp[4] + dNudy[17]*tau_gp[5] + dNudz[17]*tau_gp[2] - dNudz[17]*p_gp;
  Y[54] += dNudx[18]*tau_gp[0] + dNudy[18]*tau_gp[3] + dNudz[18]*tau_gp[4] - dNudx[18]*p_gp;
  Y[55] += dNudx[18]*tau_gp[3] + dNudy[18]*tau_gp[1] + dNudz[18]*tau_gp[5] - dNudy[18]*p_gp;
  Y[56] += dNudx[18]*tau_gp[4] + dNudy[18]*tau_gp[5] + dNudz[18]*tau_gp[2] - dNudz[18]*p_gp;
  Y[57] += dNudx[19]*tau_gp[0] + dNudy[19]*tau_gp[3] + dNudz[19]*tau_gp[4] - dNudx[19]*p_gp;
  Y[58] += dNudx[19]*tau_gp[3] + dNudy[19]*tau_gp[1] + dNudz[19]*tau_gp[5] - dNudy[19]*p_gp;
  Y[59] += dNudx[19]*tau_gp[4] + dNudy[19]*tau_gp[5] + dNudz[19]*tau_gp[2] - dNudz[19]*p_gp;
  Y[60] += dNudx[20]*tau_gp[0] + dNudy[20]*tau_gp[3] + dNudz[20]*tau_gp[4] - dNudx[20]*p_gp;
  Y[61] += dNudx[20]*tau_gp[3] + dNudy[20]*tau_gp[1] + dNudz[20]*tau_gp[5] - dNudy[20]*p_gp;
  Y[62] += dNudx[20]*tau_gp[4] + dNudy[20]*tau_gp[5] + dNudz[20]*tau_gp[2] - dNudz[20]*p_gp;
  Y[63] += dNudx[21]*tau_gp[0] + dNudy[21]*tau_gp[3] + dNudz[21]*tau_gp[4] - dNudx[21]*p_gp;
  Y[64] += dNudx[21]*tau_gp[3] + dNudy[21]*tau_gp[1] + dNudz[21]*tau_gp[5] - dNudy[21]*p_gp;
  Y[65] += dNudx[21]*tau_gp[4] + dNudy[21]*tau_gp[5] + dNudz[21]*tau_gp[2] - dNudz[21]*p_gp;
  Y[66] += dNudx[22]*tau_gp[0] + dNudy[22]*tau_gp[3] + dNudz[22]*tau_gp[4] - dNudx[22]*p_gp;
  Y[67] += dNudx[22]*tau_gp[3] + dNudy[22]*tau_gp[1] + dNudz[22]*tau_gp[5] - dNudy[22]*p_gp;
  Y[68] += dNudx[22]*tau_gp[4] + dNudy[22]*tau_gp[5] + dNudz[22]*tau_gp[2] - dNudz[22]*p_gp;
  Y[69] += dNudx[23]*tau_gp[0] + dNudy[23]*tau_gp[3] + dNudz[23]*tau_gp[4] - dNudx[23]*p_gp;
  Y[70] += dNudx[23]*tau_gp[3] + dNudy[23]*tau_gp[1] + dNudz[23]*tau_gp[5] - dNudy[23]*p_gp;
  Y[71] += dNudx[23]*tau_gp[4] + dNudy[23]*tau_gp[5] + dNudz[23]*tau_gp[2] - dNudz[23]*p_gp;
  Y[72] += dNudx[24]*tau_gp[0] + dNudy[24]*tau_gp[3] + dNudz[24]*tau_gp[4] - dNudx[24]*p_gp;
  Y[73] += dNudx[24]*tau_gp[3] + dNudy[24]*tau_gp[1] + dNudz[24]*tau_gp[5] - dNudy[24]*p_gp;
  Y[74] += dNudx[24]*tau_gp[4] + dNudy[24]*tau_gp[5] + dNudz[24]*tau_gp[2] - dNudz[24]*p_gp;
  Y[75] += dNudx[25]*tau_gp[0] + dNudy[25]*tau_gp[3] + dNudz[25]*tau_gp[4] - dNudx[25]*p_gp;
  Y[76] += dNudx[25]*tau_gp[3] + dNudy[25]*tau_gp[1] + dNudz[25]*tau_gp[5] - dNudy[25]*p_gp;
  Y[77] += dNudx[25]*tau_gp[4] + dNudy[25]*tau_gp[5] + dNudz[25]*tau_gp[2] - dNudz[25]*p_gp;
  Y[78] += dNudx[26]*tau_gp[0] + dNudy[26]*tau_gp[3] + dNudz[26]*tau_gp[4] - dNudx[26]*p_gp;
  Y[79] += dNudx[26]*tau_gp[3] + dNudy[26]*tau_gp[1] + dNudz[26]*tau_gp[5] - dNudy[26]*p_gp;
  Y[80] += dNudx[26]*tau_gp[4] + dNudy[26]*tau_gp[5] + dNudz[26]*tau_gp[2] - dNudz[26]*p_gp;
  Y[81] += -Np[0]*div_gp;
  Y[82] += -Np[1]*div_gp;
  Y[83] += -Np[2]*div_gp;
  Y[84] += -Np[3]*div_gp;
  
// total operations = 1158
}
inline void MatMultMF_Stokes_MixedFEM3d_B11(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = 0;
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = Ux[0]*dNudx[0] + Ux[10]*dNudx[10] + Ux[11]*dNudx[11] + Ux[12]*dNudx[12] + Ux[13]*dNudx[13] + Ux[14]*dNudx[14] + Ux[15]*dNudx[15] + Ux[16]*dNudx[16] + Ux[17]*dNudx[17] + Ux[18]*dNudx[18] + Ux[19]*dNudx[19] + Ux[1]*dNudx[1] + Ux[20]*dNudx[20] + Ux[21]*dNudx[21] + Ux[22]*dNudx[22] + Ux[23]*dNudx[23] + Ux[24]*dNudx[24] + Ux[25]*dNudx[25] + Ux[26]*dNudx[26] + Ux[2]*dNudx[2] + Ux[3]*dNudx[3] + Ux[4]*dNudx[4] + Ux[5]*dNudx[5] + Ux[6]*dNudx[6] + Ux[7]*dNudx[7] + Ux[8]*dNudx[8] + Ux[9]*dNudx[9];
  strain_rate_gp[1] = Uy[0]*dNudy[0] + Uy[10]*dNudy[10] + Uy[11]*dNudy[11] + Uy[12]*dNudy[12] + Uy[13]*dNudy[13] + Uy[14]*dNudy[14] + Uy[15]*dNudy[15] + Uy[16]*dNudy[16] + Uy[17]*dNudy[17] + Uy[18]*dNudy[18] + Uy[19]*dNudy[19] + Uy[1]*dNudy[1] + Uy[20]*dNudy[20] + Uy[21]*dNudy[21] + Uy[22]*dNudy[22] + Uy[23]*dNudy[23] + Uy[24]*dNudy[24] + Uy[25]*dNudy[25] + Uy[26]*dNudy[26] + Uy[2]*dNudy[2] + Uy[3]*dNudy[3] + Uy[4]*dNudy[4] + Uy[5]*dNudy[5] + Uy[6]*dNudy[6] + Uy[7]*dNudy[7] + Uy[8]*dNudy[8] + Uy[9]*dNudy[9];
  strain_rate_gp[2] = Uz[0]*dNudz[0] + Uz[10]*dNudz[10] + Uz[11]*dNudz[11] + Uz[12]*dNudz[12] + Uz[13]*dNudz[13] + Uz[14]*dNudz[14] + Uz[15]*dNudz[15] + Uz[16]*dNudz[16] + Uz[17]*dNudz[17] + Uz[18]*dNudz[18] + Uz[19]*dNudz[19] + Uz[1]*dNudz[1] + Uz[20]*dNudz[20] + Uz[21]*dNudz[21] + Uz[22]*dNudz[22] + Uz[23]*dNudz[23] + Uz[24]*dNudz[24] + Uz[25]*dNudz[25] + Uz[26]*dNudz[26] + Uz[2]*dNudz[2] + Uz[3]*dNudz[3] + Uz[4]*dNudz[4] + Uz[5]*dNudz[5] + Uz[6]*dNudz[6] + Uz[7]*dNudz[7] + Uz[8]*dNudz[8] + Uz[9]*dNudz[9];
  strain_rate_gp[3] = Ux[0]*dNudy[0] + Ux[10]*dNudy[10] + Ux[11]*dNudy[11] + Ux[12]*dNudy[12] + Ux[13]*dNudy[13] + Ux[14]*dNudy[14] + Ux[15]*dNudy[15] + Ux[16]*dNudy[16] + Ux[17]*dNudy[17] + Ux[18]*dNudy[18] + Ux[19]*dNudy[19] + Ux[1]*dNudy[1] + Ux[20]*dNudy[20] + Ux[21]*dNudy[21] + Ux[22]*dNudy[22] + Ux[23]*dNudy[23] + Ux[24]*dNudy[24] + Ux[25]*dNudy[25] + Ux[26]*dNudy[26] + Ux[2]*dNudy[2] + Ux[3]*dNudy[3] + Ux[4]*dNudy[4] + Ux[5]*dNudy[5] + Ux[6]*dNudy[6] + Ux[7]*dNudy[7] + Ux[8]*dNudy[8] + Ux[9]*dNudy[9] + Uy[0]*dNudx[0] + Uy[10]*dNudx[10] + Uy[11]*dNudx[11] + Uy[12]*dNudx[12] + Uy[13]*dNudx[13] + Uy[14]*dNudx[14] + Uy[15]*dNudx[15] + Uy[16]*dNudx[16] + Uy[17]*dNudx[17] + Uy[18]*dNudx[18] + Uy[19]*dNudx[19] + Uy[1]*dNudx[1] + Uy[20]*dNudx[20] + Uy[21]*dNudx[21] + Uy[22]*dNudx[22] + Uy[23]*dNudx[23] + Uy[24]*dNudx[24] + Uy[25]*dNudx[25] + Uy[26]*dNudx[26] + Uy[2]*dNudx[2] + Uy[3]*dNudx[3] + Uy[4]*dNudx[4] + Uy[5]*dNudx[5] + Uy[6]*dNudx[6] + Uy[7]*dNudx[7] + Uy[8]*dNudx[8] + Uy[9]*dNudx[9];
  strain_rate_gp[4] = Ux[0]*dNudz[0] + Ux[10]*dNudz[10] + Ux[11]*dNudz[11] + Ux[12]*dNudz[12] + Ux[13]*dNudz[13] + Ux[14]*dNudz[14] + Ux[15]*dNudz[15] + Ux[16]*dNudz[16] + Ux[17]*dNudz[17] + Ux[18]*dNudz[18] + Ux[19]*dNudz[19] + Ux[1]*dNudz[1] + Ux[20]*dNudz[20] + Ux[21]*dNudz[21] + Ux[22]*dNudz[22] + Ux[23]*dNudz[23] + Ux[24]*dNudz[24] + Ux[25]*dNudz[25] + Ux[26]*dNudz[26] + Ux[2]*dNudz[2] + Ux[3]*dNudz[3] + Ux[4]*dNudz[4] + Ux[5]*dNudz[5] + Ux[6]*dNudz[6] + Ux[7]*dNudz[7] + Ux[8]*dNudz[8] + Ux[9]*dNudz[9] + Uz[0]*dNudx[0] + Uz[10]*dNudx[10] + Uz[11]*dNudx[11] + Uz[12]*dNudx[12] + Uz[13]*dNudx[13] + Uz[14]*dNudx[14] + Uz[15]*dNudx[15] + Uz[16]*dNudx[16] + Uz[17]*dNudx[17] + Uz[18]*dNudx[18] + Uz[19]*dNudx[19] + Uz[1]*dNudx[1] + Uz[20]*dNudx[20] + Uz[21]*dNudx[21] + Uz[22]*dNudx[22] + Uz[23]*dNudx[23] + Uz[24]*dNudx[24] + Uz[25]*dNudx[25] + Uz[26]*dNudx[26] + Uz[2]*dNudx[2] + Uz[3]*dNudx[3] + Uz[4]*dNudx[4] + Uz[5]*dNudx[5] + Uz[6]*dNudx[6] + Uz[7]*dNudx[7] + Uz[8]*dNudx[8] + Uz[9]*dNudx[9];
  strain_rate_gp[5] = Uy[0]*dNudz[0] + Uy[10]*dNudz[10] + Uy[11]*dNudz[11] + Uy[12]*dNudz[12] + Uy[13]*dNudz[13] + Uy[14]*dNudz[14] + Uy[15]*dNudz[15] + Uy[16]*dNudz[16] + Uy[17]*dNudz[17] + Uy[18]*dNudz[18] + Uy[19]*dNudz[19] + Uy[1]*dNudz[1] + Uy[20]*dNudz[20] + Uy[21]*dNudz[21] + Uy[22]*dNudz[22] + Uy[23]*dNudz[23] + Uy[24]*dNudz[24] + Uy[25]*dNudz[25] + Uy[26]*dNudz[26] + Uy[2]*dNudz[2] + Uy[3]*dNudz[3] + Uy[4]*dNudz[4] + Uy[5]*dNudz[5] + Uy[6]*dNudz[6] + Uy[7]*dNudz[7] + Uy[8]*dNudz[8] + Uy[9]*dNudz[9] + Uz[0]*dNudy[0] + Uz[10]*dNudy[10] + Uz[11]*dNudy[11] + Uz[12]*dNudy[12] + Uz[13]*dNudy[13] + Uz[14]*dNudy[14] + Uz[15]*dNudy[15] + Uz[16]*dNudy[16] + Uz[17]*dNudy[17] + Uz[18]*dNudy[18] + Uz[19]*dNudy[19] + Uz[1]*dNudy[1] + Uz[20]*dNudy[20] + Uz[21]*dNudy[21] + Uz[22]*dNudy[22] + Uz[23]*dNudy[23] + Uz[24]*dNudy[24] + Uz[25]*dNudy[25] + Uz[26]*dNudy[26] + Uz[2]*dNudy[2] + Uz[3]*dNudy[3] + Uz[4]*dNudy[4] + Uz[5]*dNudy[5] + Uz[6]*dNudy[6] + Uz[7]*dNudy[7] + Uz[8]*dNudy[8] + Uz[9]*dNudy[9];
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 2.0*eta_gp*strain_rate_gp[0];
  tau_gp[1] = 2.0*eta_gp*strain_rate_gp[1];
  tau_gp[2] = 2.0*eta_gp*strain_rate_gp[2];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[0] + strain_rate_gp[1] + strain_rate_gp[2];
  // y = A21.u at gauss point //
  Y[0] += dNudx[0]*tau_gp[0] + dNudy[0]*tau_gp[3] + dNudz[0]*tau_gp[4];
  Y[1] += dNudx[0]*tau_gp[3] + dNudy[0]*tau_gp[1] + dNudz[0]*tau_gp[5];
  Y[2] += dNudx[0]*tau_gp[4] + dNudy[0]*tau_gp[5] + dNudz[0]*tau_gp[2];
  Y[3] += dNudx[1]*tau_gp[0] + dNudy[1]*tau_gp[3] + dNudz[1]*tau_gp[4];
  Y[4] += dNudx[1]*tau_gp[3] + dNudy[1]*tau_gp[1] + dNudz[1]*tau_gp[5];
  Y[5] += dNudx[1]*tau_gp[4] + dNudy[1]*tau_gp[5] + dNudz[1]*tau_gp[2];
  Y[6] += dNudx[2]*tau_gp[0] + dNudy[2]*tau_gp[3] + dNudz[2]*tau_gp[4];
  Y[7] += dNudx[2]*tau_gp[3] + dNudy[2]*tau_gp[1] + dNudz[2]*tau_gp[5];
  Y[8] += dNudx[2]*tau_gp[4] + dNudy[2]*tau_gp[5] + dNudz[2]*tau_gp[2];
  Y[9] += dNudx[3]*tau_gp[0] + dNudy[3]*tau_gp[3] + dNudz[3]*tau_gp[4];
  Y[10] += dNudx[3]*tau_gp[3] + dNudy[3]*tau_gp[1] + dNudz[3]*tau_gp[5];
  Y[11] += dNudx[3]*tau_gp[4] + dNudy[3]*tau_gp[5] + dNudz[3]*tau_gp[2];
  Y[12] += dNudx[4]*tau_gp[0] + dNudy[4]*tau_gp[3] + dNudz[4]*tau_gp[4];
  Y[13] += dNudx[4]*tau_gp[3] + dNudy[4]*tau_gp[1] + dNudz[4]*tau_gp[5];
  Y[14] += dNudx[4]*tau_gp[4] + dNudy[4]*tau_gp[5] + dNudz[4]*tau_gp[2];
  Y[15] += dNudx[5]*tau_gp[0] + dNudy[5]*tau_gp[3] + dNudz[5]*tau_gp[4];
  Y[16] += dNudx[5]*tau_gp[3] + dNudy[5]*tau_gp[1] + dNudz[5]*tau_gp[5];
  Y[17] += dNudx[5]*tau_gp[4] + dNudy[5]*tau_gp[5] + dNudz[5]*tau_gp[2];
  Y[18] += dNudx[6]*tau_gp[0] + dNudy[6]*tau_gp[3] + dNudz[6]*tau_gp[4];
  Y[19] += dNudx[6]*tau_gp[3] + dNudy[6]*tau_gp[1] + dNudz[6]*tau_gp[5];
  Y[20] += dNudx[6]*tau_gp[4] + dNudy[6]*tau_gp[5] + dNudz[6]*tau_gp[2];
  Y[21] += dNudx[7]*tau_gp[0] + dNudy[7]*tau_gp[3] + dNudz[7]*tau_gp[4];
  Y[22] += dNudx[7]*tau_gp[3] + dNudy[7]*tau_gp[1] + dNudz[7]*tau_gp[5];
  Y[23] += dNudx[7]*tau_gp[4] + dNudy[7]*tau_gp[5] + dNudz[7]*tau_gp[2];
  Y[24] += dNudx[8]*tau_gp[0] + dNudy[8]*tau_gp[3] + dNudz[8]*tau_gp[4];
  Y[25] += dNudx[8]*tau_gp[3] + dNudy[8]*tau_gp[1] + dNudz[8]*tau_gp[5];
  Y[26] += dNudx[8]*tau_gp[4] + dNudy[8]*tau_gp[5] + dNudz[8]*tau_gp[2];
  Y[27] += dNudx[9]*tau_gp[0] + dNudy[9]*tau_gp[3] + dNudz[9]*tau_gp[4];
  Y[28] += dNudx[9]*tau_gp[3] + dNudy[9]*tau_gp[1] + dNudz[9]*tau_gp[5];
  Y[29] += dNudx[9]*tau_gp[4] + dNudy[9]*tau_gp[5] + dNudz[9]*tau_gp[2];
  Y[30] += dNudx[10]*tau_gp[0] + dNudy[10]*tau_gp[3] + dNudz[10]*tau_gp[4];
  Y[31] += dNudx[10]*tau_gp[3] + dNudy[10]*tau_gp[1] + dNudz[10]*tau_gp[5];
  Y[32] += dNudx[10]*tau_gp[4] + dNudy[10]*tau_gp[5] + dNudz[10]*tau_gp[2];
  Y[33] += dNudx[11]*tau_gp[0] + dNudy[11]*tau_gp[3] + dNudz[11]*tau_gp[4];
  Y[34] += dNudx[11]*tau_gp[3] + dNudy[11]*tau_gp[1] + dNudz[11]*tau_gp[5];
  Y[35] += dNudx[11]*tau_gp[4] + dNudy[11]*tau_gp[5] + dNudz[11]*tau_gp[2];
  Y[36] += dNudx[12]*tau_gp[0] + dNudy[12]*tau_gp[3] + dNudz[12]*tau_gp[4];
  Y[37] += dNudx[12]*tau_gp[3] + dNudy[12]*tau_gp[1] + dNudz[12]*tau_gp[5];
  Y[38] += dNudx[12]*tau_gp[4] + dNudy[12]*tau_gp[5] + dNudz[12]*tau_gp[2];
  Y[39] += dNudx[13]*tau_gp[0] + dNudy[13]*tau_gp[3] + dNudz[13]*tau_gp[4];
  Y[40] += dNudx[13]*tau_gp[3] + dNudy[13]*tau_gp[1] + dNudz[13]*tau_gp[5];
  Y[41] += dNudx[13]*tau_gp[4] + dNudy[13]*tau_gp[5] + dNudz[13]*tau_gp[2];
  Y[42] += dNudx[14]*tau_gp[0] + dNudy[14]*tau_gp[3] + dNudz[14]*tau_gp[4];
  Y[43] += dNudx[14]*tau_gp[3] + dNudy[14]*tau_gp[1] + dNudz[14]*tau_gp[5];
  Y[44] += dNudx[14]*tau_gp[4] + dNudy[14]*tau_gp[5] + dNudz[14]*tau_gp[2];
  Y[45] += dNudx[15]*tau_gp[0] + dNudy[15]*tau_gp[3] + dNudz[15]*tau_gp[4];
  Y[46] += dNudx[15]*tau_gp[3] + dNudy[15]*tau_gp[1] + dNudz[15]*tau_gp[5];
  Y[47] += dNudx[15]*tau_gp[4] + dNudy[15]*tau_gp[5] + dNudz[15]*tau_gp[2];
  Y[48] += dNudx[16]*tau_gp[0] + dNudy[16]*tau_gp[3] + dNudz[16]*tau_gp[4];
  Y[49] += dNudx[16]*tau_gp[3] + dNudy[16]*tau_gp[1] + dNudz[16]*tau_gp[5];
  Y[50] += dNudx[16]*tau_gp[4] + dNudy[16]*tau_gp[5] + dNudz[16]*tau_gp[2];
  Y[51] += dNudx[17]*tau_gp[0] + dNudy[17]*tau_gp[3] + dNudz[17]*tau_gp[4];
  Y[52] += dNudx[17]*tau_gp[3] + dNudy[17]*tau_gp[1] + dNudz[17]*tau_gp[5];
  Y[53] += dNudx[17]*tau_gp[4] + dNudy[17]*tau_gp[5] + dNudz[17]*tau_gp[2];
  Y[54] += dNudx[18]*tau_gp[0] + dNudy[18]*tau_gp[3] + dNudz[18]*tau_gp[4];
  Y[55] += dNudx[18]*tau_gp[3] + dNudy[18]*tau_gp[1] + dNudz[18]*tau_gp[5];
  Y[56] += dNudx[18]*tau_gp[4] + dNudy[18]*tau_gp[5] + dNudz[18]*tau_gp[2];
  Y[57] += dNudx[19]*tau_gp[0] + dNudy[19]*tau_gp[3] + dNudz[19]*tau_gp[4];
  Y[58] += dNudx[19]*tau_gp[3] + dNudy[19]*tau_gp[1] + dNudz[19]*tau_gp[5];
  Y[59] += dNudx[19]*tau_gp[4] + dNudy[19]*tau_gp[5] + dNudz[19]*tau_gp[2];
  Y[60] += dNudx[20]*tau_gp[0] + dNudy[20]*tau_gp[3] + dNudz[20]*tau_gp[4];
  Y[61] += dNudx[20]*tau_gp[3] + dNudy[20]*tau_gp[1] + dNudz[20]*tau_gp[5];
  Y[62] += dNudx[20]*tau_gp[4] + dNudy[20]*tau_gp[5] + dNudz[20]*tau_gp[2];
  Y[63] += dNudx[21]*tau_gp[0] + dNudy[21]*tau_gp[3] + dNudz[21]*tau_gp[4];
  Y[64] += dNudx[21]*tau_gp[3] + dNudy[21]*tau_gp[1] + dNudz[21]*tau_gp[5];
  Y[65] += dNudx[21]*tau_gp[4] + dNudy[21]*tau_gp[5] + dNudz[21]*tau_gp[2];
  Y[66] += dNudx[22]*tau_gp[0] + dNudy[22]*tau_gp[3] + dNudz[22]*tau_gp[4];
  Y[67] += dNudx[22]*tau_gp[3] + dNudy[22]*tau_gp[1] + dNudz[22]*tau_gp[5];
  Y[68] += dNudx[22]*tau_gp[4] + dNudy[22]*tau_gp[5] + dNudz[22]*tau_gp[2];
  Y[69] += dNudx[23]*tau_gp[0] + dNudy[23]*tau_gp[3] + dNudz[23]*tau_gp[4];
  Y[70] += dNudx[23]*tau_gp[3] + dNudy[23]*tau_gp[1] + dNudz[23]*tau_gp[5];
  Y[71] += dNudx[23]*tau_gp[4] + dNudy[23]*tau_gp[5] + dNudz[23]*tau_gp[2];
  Y[72] += dNudx[24]*tau_gp[0] + dNudy[24]*tau_gp[3] + dNudz[24]*tau_gp[4];
  Y[73] += dNudx[24]*tau_gp[3] + dNudy[24]*tau_gp[1] + dNudz[24]*tau_gp[5];
  Y[74] += dNudx[24]*tau_gp[4] + dNudy[24]*tau_gp[5] + dNudz[24]*tau_gp[2];
  Y[75] += dNudx[25]*tau_gp[0] + dNudy[25]*tau_gp[3] + dNudz[25]*tau_gp[4];
  Y[76] += dNudx[25]*tau_gp[3] + dNudy[25]*tau_gp[1] + dNudz[25]*tau_gp[5];
  Y[77] += dNudx[25]*tau_gp[4] + dNudy[25]*tau_gp[5] + dNudz[25]*tau_gp[2];
  Y[78] += dNudx[26]*tau_gp[0] + dNudy[26]*tau_gp[3] + dNudz[26]*tau_gp[4];
  Y[79] += dNudx[26]*tau_gp[3] + dNudy[26]*tau_gp[1] + dNudz[26]*tau_gp[5];
  Y[80] += dNudx[26]*tau_gp[4] + dNudy[26]*tau_gp[5] + dNudz[26]*tau_gp[2];
  
// total operations = 977
}
inline void MatMultMF_Stokes_MixedFEM3d_Buu(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = 0;
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = Ux[0]*dNudx[0] + Ux[10]*dNudx[10] + Ux[11]*dNudx[11] + Ux[12]*dNudx[12] + Ux[13]*dNudx[13] + Ux[14]*dNudx[14] + Ux[15]*dNudx[15] + Ux[16]*dNudx[16] + Ux[17]*dNudx[17] + Ux[18]*dNudx[18] + Ux[19]*dNudx[19] + Ux[1]*dNudx[1] + Ux[20]*dNudx[20] + Ux[21]*dNudx[21] + Ux[22]*dNudx[22] + Ux[23]*dNudx[23] + Ux[24]*dNudx[24] + Ux[25]*dNudx[25] + Ux[26]*dNudx[26] + Ux[2]*dNudx[2] + Ux[3]*dNudx[3] + Ux[4]*dNudx[4] + Ux[5]*dNudx[5] + Ux[6]*dNudx[6] + Ux[7]*dNudx[7] + Ux[8]*dNudx[8] + Ux[9]*dNudx[9];
  strain_rate_gp[1] = 0;
  strain_rate_gp[2] = 0;
  strain_rate_gp[3] = Ux[0]*dNudy[0] + Ux[10]*dNudy[10] + Ux[11]*dNudy[11] + Ux[12]*dNudy[12] + Ux[13]*dNudy[13] + Ux[14]*dNudy[14] + Ux[15]*dNudy[15] + Ux[16]*dNudy[16] + Ux[17]*dNudy[17] + Ux[18]*dNudy[18] + Ux[19]*dNudy[19] + Ux[1]*dNudy[1] + Ux[20]*dNudy[20] + Ux[21]*dNudy[21] + Ux[22]*dNudy[22] + Ux[23]*dNudy[23] + Ux[24]*dNudy[24] + Ux[25]*dNudy[25] + Ux[26]*dNudy[26] + Ux[2]*dNudy[2] + Ux[3]*dNudy[3] + Ux[4]*dNudy[4] + Ux[5]*dNudy[5] + Ux[6]*dNudy[6] + Ux[7]*dNudy[7] + Ux[8]*dNudy[8] + Ux[9]*dNudy[9];
  strain_rate_gp[4] = Ux[0]*dNudz[0] + Ux[10]*dNudz[10] + Ux[11]*dNudz[11] + Ux[12]*dNudz[12] + Ux[13]*dNudz[13] + Ux[14]*dNudz[14] + Ux[15]*dNudz[15] + Ux[16]*dNudz[16] + Ux[17]*dNudz[17] + Ux[18]*dNudz[18] + Ux[19]*dNudz[19] + Ux[1]*dNudz[1] + Ux[20]*dNudz[20] + Ux[21]*dNudz[21] + Ux[22]*dNudz[22] + Ux[23]*dNudz[23] + Ux[24]*dNudz[24] + Ux[25]*dNudz[25] + Ux[26]*dNudz[26] + Ux[2]*dNudz[2] + Ux[3]*dNudz[3] + Ux[4]*dNudz[4] + Ux[5]*dNudz[5] + Ux[6]*dNudz[6] + Ux[7]*dNudz[7] + Ux[8]*dNudz[8] + Ux[9]*dNudz[9];
  strain_rate_gp[5] = 0;
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 2.0*eta_gp*strain_rate_gp[0];
  tau_gp[1] = 0;
  tau_gp[2] = 0;
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = 0;
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[0];
  // y = A21.u at gauss point //
  Y[0] += dNudx[0]*tau_gp[0] + dNudy[0]*tau_gp[3] + dNudz[0]*tau_gp[4];
  Y[1] += dNudx[1]*tau_gp[0] + dNudy[1]*tau_gp[3] + dNudz[1]*tau_gp[4];
  Y[2] += dNudx[2]*tau_gp[0] + dNudy[2]*tau_gp[3] + dNudz[2]*tau_gp[4];
  Y[3] += dNudx[3]*tau_gp[0] + dNudy[3]*tau_gp[3] + dNudz[3]*tau_gp[4];
  Y[4] += dNudx[4]*tau_gp[0] + dNudy[4]*tau_gp[3] + dNudz[4]*tau_gp[4];
  Y[5] += dNudx[5]*tau_gp[0] + dNudy[5]*tau_gp[3] + dNudz[5]*tau_gp[4];
  Y[6] += dNudx[6]*tau_gp[0] + dNudy[6]*tau_gp[3] + dNudz[6]*tau_gp[4];
  Y[7] += dNudx[7]*tau_gp[0] + dNudy[7]*tau_gp[3] + dNudz[7]*tau_gp[4];
  Y[8] += dNudx[8]*tau_gp[0] + dNudy[8]*tau_gp[3] + dNudz[8]*tau_gp[4];
  Y[9] += dNudx[9]*tau_gp[0] + dNudy[9]*tau_gp[3] + dNudz[9]*tau_gp[4];
  Y[10] += dNudx[10]*tau_gp[0] + dNudy[10]*tau_gp[3] + dNudz[10]*tau_gp[4];
  Y[11] += dNudx[11]*tau_gp[0] + dNudy[11]*tau_gp[3] + dNudz[11]*tau_gp[4];
  Y[12] += dNudx[12]*tau_gp[0] + dNudy[12]*tau_gp[3] + dNudz[12]*tau_gp[4];
  Y[13] += dNudx[13]*tau_gp[0] + dNudy[13]*tau_gp[3] + dNudz[13]*tau_gp[4];
  Y[14] += dNudx[14]*tau_gp[0] + dNudy[14]*tau_gp[3] + dNudz[14]*tau_gp[4];
  Y[15] += dNudx[15]*tau_gp[0] + dNudy[15]*tau_gp[3] + dNudz[15]*tau_gp[4];
  Y[16] += dNudx[16]*tau_gp[0] + dNudy[16]*tau_gp[3] + dNudz[16]*tau_gp[4];
  Y[17] += dNudx[17]*tau_gp[0] + dNudy[17]*tau_gp[3] + dNudz[17]*tau_gp[4];
  Y[18] += dNudx[18]*tau_gp[0] + dNudy[18]*tau_gp[3] + dNudz[18]*tau_gp[4];
  Y[19] += dNudx[19]*tau_gp[0] + dNudy[19]*tau_gp[3] + dNudz[19]*tau_gp[4];
  Y[20] += dNudx[20]*tau_gp[0] + dNudy[20]*tau_gp[3] + dNudz[20]*tau_gp[4];
  Y[21] += dNudx[21]*tau_gp[0] + dNudy[21]*tau_gp[3] + dNudz[21]*tau_gp[4];
  Y[22] += dNudx[22]*tau_gp[0] + dNudy[22]*tau_gp[3] + dNudz[22]*tau_gp[4];
  Y[23] += dNudx[23]*tau_gp[0] + dNudy[23]*tau_gp[3] + dNudz[23]*tau_gp[4];
  Y[24] += dNudx[24]*tau_gp[0] + dNudy[24]*tau_gp[3] + dNudz[24]*tau_gp[4];
  Y[25] += dNudx[25]*tau_gp[0] + dNudy[25]*tau_gp[3] + dNudz[25]*tau_gp[4];
  Y[26] += dNudx[26]*tau_gp[0] + dNudy[26]*tau_gp[3] + dNudz[26]*tau_gp[4];
  
// total operations = 325
}
inline void MatMultMF_Stokes_MixedFEM3d_Bvv(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = 0;
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = 0;
  strain_rate_gp[1] = Uy[0]*dNudy[0] + Uy[10]*dNudy[10] + Uy[11]*dNudy[11] + Uy[12]*dNudy[12] + Uy[13]*dNudy[13] + Uy[14]*dNudy[14] + Uy[15]*dNudy[15] + Uy[16]*dNudy[16] + Uy[17]*dNudy[17] + Uy[18]*dNudy[18] + Uy[19]*dNudy[19] + Uy[1]*dNudy[1] + Uy[20]*dNudy[20] + Uy[21]*dNudy[21] + Uy[22]*dNudy[22] + Uy[23]*dNudy[23] + Uy[24]*dNudy[24] + Uy[25]*dNudy[25] + Uy[26]*dNudy[26] + Uy[2]*dNudy[2] + Uy[3]*dNudy[3] + Uy[4]*dNudy[4] + Uy[5]*dNudy[5] + Uy[6]*dNudy[6] + Uy[7]*dNudy[7] + Uy[8]*dNudy[8] + Uy[9]*dNudy[9];
  strain_rate_gp[2] = 0;
  strain_rate_gp[3] = Uy[0]*dNudx[0] + Uy[10]*dNudx[10] + Uy[11]*dNudx[11] + Uy[12]*dNudx[12] + Uy[13]*dNudx[13] + Uy[14]*dNudx[14] + Uy[15]*dNudx[15] + Uy[16]*dNudx[16] + Uy[17]*dNudx[17] + Uy[18]*dNudx[18] + Uy[19]*dNudx[19] + Uy[1]*dNudx[1] + Uy[20]*dNudx[20] + Uy[21]*dNudx[21] + Uy[22]*dNudx[22] + Uy[23]*dNudx[23] + Uy[24]*dNudx[24] + Uy[25]*dNudx[25] + Uy[26]*dNudx[26] + Uy[2]*dNudx[2] + Uy[3]*dNudx[3] + Uy[4]*dNudx[4] + Uy[5]*dNudx[5] + Uy[6]*dNudx[6] + Uy[7]*dNudx[7] + Uy[8]*dNudx[8] + Uy[9]*dNudx[9];
  strain_rate_gp[4] = 0;
  strain_rate_gp[5] = Uy[0]*dNudz[0] + Uy[10]*dNudz[10] + Uy[11]*dNudz[11] + Uy[12]*dNudz[12] + Uy[13]*dNudz[13] + Uy[14]*dNudz[14] + Uy[15]*dNudz[15] + Uy[16]*dNudz[16] + Uy[17]*dNudz[17] + Uy[18]*dNudz[18] + Uy[19]*dNudz[19] + Uy[1]*dNudz[1] + Uy[20]*dNudz[20] + Uy[21]*dNudz[21] + Uy[22]*dNudz[22] + Uy[23]*dNudz[23] + Uy[24]*dNudz[24] + Uy[25]*dNudz[25] + Uy[26]*dNudz[26] + Uy[2]*dNudz[2] + Uy[3]*dNudz[3] + Uy[4]*dNudz[4] + Uy[5]*dNudz[5] + Uy[6]*dNudz[6] + Uy[7]*dNudz[7] + Uy[8]*dNudz[8] + Uy[9]*dNudz[9];
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 0;
  tau_gp[1] = 2.0*eta_gp*strain_rate_gp[1];
  tau_gp[2] = 0;
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = 0;
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[1];
  // y = A21.u at gauss point //
  Y[0] += dNudx[0]*tau_gp[3] + dNudy[0]*tau_gp[1] + dNudz[0]*tau_gp[5];
  Y[1] += dNudx[1]*tau_gp[3] + dNudy[1]*tau_gp[1] + dNudz[1]*tau_gp[5];
  Y[2] += dNudx[2]*tau_gp[3] + dNudy[2]*tau_gp[1] + dNudz[2]*tau_gp[5];
  Y[3] += dNudx[3]*tau_gp[3] + dNudy[3]*tau_gp[1] + dNudz[3]*tau_gp[5];
  Y[4] += dNudx[4]*tau_gp[3] + dNudy[4]*tau_gp[1] + dNudz[4]*tau_gp[5];
  Y[5] += dNudx[5]*tau_gp[3] + dNudy[5]*tau_gp[1] + dNudz[5]*tau_gp[5];
  Y[6] += dNudx[6]*tau_gp[3] + dNudy[6]*tau_gp[1] + dNudz[6]*tau_gp[5];
  Y[7] += dNudx[7]*tau_gp[3] + dNudy[7]*tau_gp[1] + dNudz[7]*tau_gp[5];
  Y[8] += dNudx[8]*tau_gp[3] + dNudy[8]*tau_gp[1] + dNudz[8]*tau_gp[5];
  Y[9] += dNudx[9]*tau_gp[3] + dNudy[9]*tau_gp[1] + dNudz[9]*tau_gp[5];
  Y[10] += dNudx[10]*tau_gp[3] + dNudy[10]*tau_gp[1] + dNudz[10]*tau_gp[5];
  Y[11] += dNudx[11]*tau_gp[3] + dNudy[11]*tau_gp[1] + dNudz[11]*tau_gp[5];
  Y[12] += dNudx[12]*tau_gp[3] + dNudy[12]*tau_gp[1] + dNudz[12]*tau_gp[5];
  Y[13] += dNudx[13]*tau_gp[3] + dNudy[13]*tau_gp[1] + dNudz[13]*tau_gp[5];
  Y[14] += dNudx[14]*tau_gp[3] + dNudy[14]*tau_gp[1] + dNudz[14]*tau_gp[5];
  Y[15] += dNudx[15]*tau_gp[3] + dNudy[15]*tau_gp[1] + dNudz[15]*tau_gp[5];
  Y[16] += dNudx[16]*tau_gp[3] + dNudy[16]*tau_gp[1] + dNudz[16]*tau_gp[5];
  Y[17] += dNudx[17]*tau_gp[3] + dNudy[17]*tau_gp[1] + dNudz[17]*tau_gp[5];
  Y[18] += dNudx[18]*tau_gp[3] + dNudy[18]*tau_gp[1] + dNudz[18]*tau_gp[5];
  Y[19] += dNudx[19]*tau_gp[3] + dNudy[19]*tau_gp[1] + dNudz[19]*tau_gp[5];
  Y[20] += dNudx[20]*tau_gp[3] + dNudy[20]*tau_gp[1] + dNudz[20]*tau_gp[5];
  Y[21] += dNudx[21]*tau_gp[3] + dNudy[21]*tau_gp[1] + dNudz[21]*tau_gp[5];
  Y[22] += dNudx[22]*tau_gp[3] + dNudy[22]*tau_gp[1] + dNudz[22]*tau_gp[5];
  Y[23] += dNudx[23]*tau_gp[3] + dNudy[23]*tau_gp[1] + dNudz[23]*tau_gp[5];
  Y[24] += dNudx[24]*tau_gp[3] + dNudy[24]*tau_gp[1] + dNudz[24]*tau_gp[5];
  Y[25] += dNudx[25]*tau_gp[3] + dNudy[25]*tau_gp[1] + dNudz[25]*tau_gp[5];
  Y[26] += dNudx[26]*tau_gp[3] + dNudy[26]*tau_gp[1] + dNudz[26]*tau_gp[5];
  
// total operations = 325
}
inline void MatMultMF_Stokes_MixedFEM3d_Bww(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = 0;
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = 0;
  strain_rate_gp[1] = 0;
  strain_rate_gp[2] = Uz[0]*dNudz[0] + Uz[10]*dNudz[10] + Uz[11]*dNudz[11] + Uz[12]*dNudz[12] + Uz[13]*dNudz[13] + Uz[14]*dNudz[14] + Uz[15]*dNudz[15] + Uz[16]*dNudz[16] + Uz[17]*dNudz[17] + Uz[18]*dNudz[18] + Uz[19]*dNudz[19] + Uz[1]*dNudz[1] + Uz[20]*dNudz[20] + Uz[21]*dNudz[21] + Uz[22]*dNudz[22] + Uz[23]*dNudz[23] + Uz[24]*dNudz[24] + Uz[25]*dNudz[25] + Uz[26]*dNudz[26] + Uz[2]*dNudz[2] + Uz[3]*dNudz[3] + Uz[4]*dNudz[4] + Uz[5]*dNudz[5] + Uz[6]*dNudz[6] + Uz[7]*dNudz[7] + Uz[8]*dNudz[8] + Uz[9]*dNudz[9];
  strain_rate_gp[3] = 0;
  strain_rate_gp[4] = Uz[0]*dNudx[0] + Uz[10]*dNudx[10] + Uz[11]*dNudx[11] + Uz[12]*dNudx[12] + Uz[13]*dNudx[13] + Uz[14]*dNudx[14] + Uz[15]*dNudx[15] + Uz[16]*dNudx[16] + Uz[17]*dNudx[17] + Uz[18]*dNudx[18] + Uz[19]*dNudx[19] + Uz[1]*dNudx[1] + Uz[20]*dNudx[20] + Uz[21]*dNudx[21] + Uz[22]*dNudx[22] + Uz[23]*dNudx[23] + Uz[24]*dNudx[24] + Uz[25]*dNudx[25] + Uz[26]*dNudx[26] + Uz[2]*dNudx[2] + Uz[3]*dNudx[3] + Uz[4]*dNudx[4] + Uz[5]*dNudx[5] + Uz[6]*dNudx[6] + Uz[7]*dNudx[7] + Uz[8]*dNudx[8] + Uz[9]*dNudx[9];
  strain_rate_gp[5] = Uz[0]*dNudy[0] + Uz[10]*dNudy[10] + Uz[11]*dNudy[11] + Uz[12]*dNudy[12] + Uz[13]*dNudy[13] + Uz[14]*dNudy[14] + Uz[15]*dNudy[15] + Uz[16]*dNudy[16] + Uz[17]*dNudy[17] + Uz[18]*dNudy[18] + Uz[19]*dNudy[19] + Uz[1]*dNudy[1] + Uz[20]*dNudy[20] + Uz[21]*dNudy[21] + Uz[22]*dNudy[22] + Uz[23]*dNudy[23] + Uz[24]*dNudy[24] + Uz[25]*dNudy[25] + Uz[26]*dNudy[26] + Uz[2]*dNudy[2] + Uz[3]*dNudy[3] + Uz[4]*dNudy[4] + Uz[5]*dNudy[5] + Uz[6]*dNudy[6] + Uz[7]*dNudy[7] + Uz[8]*dNudy[8] + Uz[9]*dNudy[9];
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 0;
  tau_gp[1] = 0;
  tau_gp[2] = 2.0*eta_gp*strain_rate_gp[2];
  tau_gp[3] = 0;
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[2];
  // y = A21.u at gauss point //
  Y[0] += dNudx[0]*tau_gp[4] + dNudy[0]*tau_gp[5] + dNudz[0]*tau_gp[2];
  Y[1] += dNudx[1]*tau_gp[4] + dNudy[1]*tau_gp[5] + dNudz[1]*tau_gp[2];
  Y[2] += dNudx[2]*tau_gp[4] + dNudy[2]*tau_gp[5] + dNudz[2]*tau_gp[2];
  Y[3] += dNudx[3]*tau_gp[4] + dNudy[3]*tau_gp[5] + dNudz[3]*tau_gp[2];
  Y[4] += dNudx[4]*tau_gp[4] + dNudy[4]*tau_gp[5] + dNudz[4]*tau_gp[2];
  Y[5] += dNudx[5]*tau_gp[4] + dNudy[5]*tau_gp[5] + dNudz[5]*tau_gp[2];
  Y[6] += dNudx[6]*tau_gp[4] + dNudy[6]*tau_gp[5] + dNudz[6]*tau_gp[2];
  Y[7] += dNudx[7]*tau_gp[4] + dNudy[7]*tau_gp[5] + dNudz[7]*tau_gp[2];
  Y[8] += dNudx[8]*tau_gp[4] + dNudy[8]*tau_gp[5] + dNudz[8]*tau_gp[2];
  Y[9] += dNudx[9]*tau_gp[4] + dNudy[9]*tau_gp[5] + dNudz[9]*tau_gp[2];
  Y[10] += dNudx[10]*tau_gp[4] + dNudy[10]*tau_gp[5] + dNudz[10]*tau_gp[2];
  Y[11] += dNudx[11]*tau_gp[4] + dNudy[11]*tau_gp[5] + dNudz[11]*tau_gp[2];
  Y[12] += dNudx[12]*tau_gp[4] + dNudy[12]*tau_gp[5] + dNudz[12]*tau_gp[2];
  Y[13] += dNudx[13]*tau_gp[4] + dNudy[13]*tau_gp[5] + dNudz[13]*tau_gp[2];
  Y[14] += dNudx[14]*tau_gp[4] + dNudy[14]*tau_gp[5] + dNudz[14]*tau_gp[2];
  Y[15] += dNudx[15]*tau_gp[4] + dNudy[15]*tau_gp[5] + dNudz[15]*tau_gp[2];
  Y[16] += dNudx[16]*tau_gp[4] + dNudy[16]*tau_gp[5] + dNudz[16]*tau_gp[2];
  Y[17] += dNudx[17]*tau_gp[4] + dNudy[17]*tau_gp[5] + dNudz[17]*tau_gp[2];
  Y[18] += dNudx[18]*tau_gp[4] + dNudy[18]*tau_gp[5] + dNudz[18]*tau_gp[2];
  Y[19] += dNudx[19]*tau_gp[4] + dNudy[19]*tau_gp[5] + dNudz[19]*tau_gp[2];
  Y[20] += dNudx[20]*tau_gp[4] + dNudy[20]*tau_gp[5] + dNudz[20]*tau_gp[2];
  Y[21] += dNudx[21]*tau_gp[4] + dNudy[21]*tau_gp[5] + dNudz[21]*tau_gp[2];
  Y[22] += dNudx[22]*tau_gp[4] + dNudy[22]*tau_gp[5] + dNudz[22]*tau_gp[2];
  Y[23] += dNudx[23]*tau_gp[4] + dNudy[23]*tau_gp[5] + dNudz[23]*tau_gp[2];
  Y[24] += dNudx[24]*tau_gp[4] + dNudy[24]*tau_gp[5] + dNudz[24]*tau_gp[2];
  Y[25] += dNudx[25]*tau_gp[4] + dNudy[25]*tau_gp[5] + dNudz[25]*tau_gp[2];
  Y[26] += dNudx[26]*tau_gp[4] + dNudy[26]*tau_gp[5] + dNudz[26]*tau_gp[2];
  
// total operations = 325
}
inline void MatMultMF_Stokes_MixedFEM3d_A12(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = Np[0]*P[0] + Np[1]*P[1] + Np[2]*P[2] + Np[3]*P[3];
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = 0;
  strain_rate_gp[1] = 0;
  strain_rate_gp[2] = 0;
  strain_rate_gp[3] = 0;
  strain_rate_gp[4] = 0;
  strain_rate_gp[5] = 0;
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 0;
  tau_gp[1] = 0;
  tau_gp[2] = 0;
  tau_gp[3] = 0;
  tau_gp[4] = 0;
  tau_gp[5] = 0;
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = 0;
  // y = A21.u at gauss point //
  Y[0] += -dNudx[0]*p_gp;
  Y[1] += -dNudy[0]*p_gp;
  Y[2] += -dNudz[0]*p_gp;
  Y[3] += -dNudx[1]*p_gp;
  Y[4] += -dNudy[1]*p_gp;
  Y[5] += -dNudz[1]*p_gp;
  Y[6] += -dNudx[2]*p_gp;
  Y[7] += -dNudy[2]*p_gp;
  Y[8] += -dNudz[2]*p_gp;
  Y[9] += -dNudx[3]*p_gp;
  Y[10] += -dNudy[3]*p_gp;
  Y[11] += -dNudz[3]*p_gp;
  Y[12] += -dNudx[4]*p_gp;
  Y[13] += -dNudy[4]*p_gp;
  Y[14] += -dNudz[4]*p_gp;
  Y[15] += -dNudx[5]*p_gp;
  Y[16] += -dNudy[5]*p_gp;
  Y[17] += -dNudz[5]*p_gp;
  Y[18] += -dNudx[6]*p_gp;
  Y[19] += -dNudy[6]*p_gp;
  Y[20] += -dNudz[6]*p_gp;
  Y[21] += -dNudx[7]*p_gp;
  Y[22] += -dNudy[7]*p_gp;
  Y[23] += -dNudz[7]*p_gp;
  Y[24] += -dNudx[8]*p_gp;
  Y[25] += -dNudy[8]*p_gp;
  Y[26] += -dNudz[8]*p_gp;
  Y[27] += -dNudx[9]*p_gp;
  Y[28] += -dNudy[9]*p_gp;
  Y[29] += -dNudz[9]*p_gp;
  Y[30] += -dNudx[10]*p_gp;
  Y[31] += -dNudy[10]*p_gp;
  Y[32] += -dNudz[10]*p_gp;
  Y[33] += -dNudx[11]*p_gp;
  Y[34] += -dNudy[11]*p_gp;
  Y[35] += -dNudz[11]*p_gp;
  Y[36] += -dNudx[12]*p_gp;
  Y[37] += -dNudy[12]*p_gp;
  Y[38] += -dNudz[12]*p_gp;
  Y[39] += -dNudx[13]*p_gp;
  Y[40] += -dNudy[13]*p_gp;
  Y[41] += -dNudz[13]*p_gp;
  Y[42] += -dNudx[14]*p_gp;
  Y[43] += -dNudy[14]*p_gp;
  Y[44] += -dNudz[14]*p_gp;
  Y[45] += -dNudx[15]*p_gp;
  Y[46] += -dNudy[15]*p_gp;
  Y[47] += -dNudz[15]*p_gp;
  Y[48] += -dNudx[16]*p_gp;
  Y[49] += -dNudy[16]*p_gp;
  Y[50] += -dNudz[16]*p_gp;
  Y[51] += -dNudx[17]*p_gp;
  Y[52] += -dNudy[17]*p_gp;
  Y[53] += -dNudz[17]*p_gp;
  Y[54] += -dNudx[18]*p_gp;
  Y[55] += -dNudy[18]*p_gp;
  Y[56] += -dNudz[18]*p_gp;
  Y[57] += -dNudx[19]*p_gp;
  Y[58] += -dNudy[19]*p_gp;
  Y[59] += -dNudz[19]*p_gp;
  Y[60] += -dNudx[20]*p_gp;
  Y[61] += -dNudy[20]*p_gp;
  Y[62] += -dNudz[20]*p_gp;
  Y[63] += -dNudx[21]*p_gp;
  Y[64] += -dNudy[21]*p_gp;
  Y[65] += -dNudz[21]*p_gp;
  Y[66] += -dNudx[22]*p_gp;
  Y[67] += -dNudy[22]*p_gp;
  Y[68] += -dNudz[22]*p_gp;
  Y[69] += -dNudx[23]*p_gp;
  Y[70] += -dNudy[23]*p_gp;
  Y[71] += -dNudz[23]*p_gp;
  Y[72] += -dNudx[24]*p_gp;
  Y[73] += -dNudy[24]*p_gp;
  Y[74] += -dNudz[24]*p_gp;
  Y[75] += -dNudx[25]*p_gp;
  Y[76] += -dNudy[25]*p_gp;
  Y[77] += -dNudz[25]*p_gp;
  Y[78] += -dNudx[26]*p_gp;
  Y[79] += -dNudy[26]*p_gp;
  Y[80] += -dNudz[26]*p_gp;
  
// total operations = 250
}
inline void MatMultMF_Stokes_MixedFEM3d_A21(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = 0;
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = Ux[0]*dNudx[0] + Ux[10]*dNudx[10] + Ux[11]*dNudx[11] + Ux[12]*dNudx[12] + Ux[13]*dNudx[13] + Ux[14]*dNudx[14] + Ux[15]*dNudx[15] + Ux[16]*dNudx[16] + Ux[17]*dNudx[17] + Ux[18]*dNudx[18] + Ux[19]*dNudx[19] + Ux[1]*dNudx[1] + Ux[20]*dNudx[20] + Ux[21]*dNudx[21] + Ux[22]*dNudx[22] + Ux[23]*dNudx[23] + Ux[24]*dNudx[24] + Ux[25]*dNudx[25] + Ux[26]*dNudx[26] + Ux[2]*dNudx[2] + Ux[3]*dNudx[3] + Ux[4]*dNudx[4] + Ux[5]*dNudx[5] + Ux[6]*dNudx[6] + Ux[7]*dNudx[7] + Ux[8]*dNudx[8] + Ux[9]*dNudx[9];
  strain_rate_gp[1] = Uy[0]*dNudy[0] + Uy[10]*dNudy[10] + Uy[11]*dNudy[11] + Uy[12]*dNudy[12] + Uy[13]*dNudy[13] + Uy[14]*dNudy[14] + Uy[15]*dNudy[15] + Uy[16]*dNudy[16] + Uy[17]*dNudy[17] + Uy[18]*dNudy[18] + Uy[19]*dNudy[19] + Uy[1]*dNudy[1] + Uy[20]*dNudy[20] + Uy[21]*dNudy[21] + Uy[22]*dNudy[22] + Uy[23]*dNudy[23] + Uy[24]*dNudy[24] + Uy[25]*dNudy[25] + Uy[26]*dNudy[26] + Uy[2]*dNudy[2] + Uy[3]*dNudy[3] + Uy[4]*dNudy[4] + Uy[5]*dNudy[5] + Uy[6]*dNudy[6] + Uy[7]*dNudy[7] + Uy[8]*dNudy[8] + Uy[9]*dNudy[9];
  strain_rate_gp[2] = Uz[0]*dNudz[0] + Uz[10]*dNudz[10] + Uz[11]*dNudz[11] + Uz[12]*dNudz[12] + Uz[13]*dNudz[13] + Uz[14]*dNudz[14] + Uz[15]*dNudz[15] + Uz[16]*dNudz[16] + Uz[17]*dNudz[17] + Uz[18]*dNudz[18] + Uz[19]*dNudz[19] + Uz[1]*dNudz[1] + Uz[20]*dNudz[20] + Uz[21]*dNudz[21] + Uz[22]*dNudz[22] + Uz[23]*dNudz[23] + Uz[24]*dNudz[24] + Uz[25]*dNudz[25] + Uz[26]*dNudz[26] + Uz[2]*dNudz[2] + Uz[3]*dNudz[3] + Uz[4]*dNudz[4] + Uz[5]*dNudz[5] + Uz[6]*dNudz[6] + Uz[7]*dNudz[7] + Uz[8]*dNudz[8] + Uz[9]*dNudz[9];
  strain_rate_gp[3] = Ux[0]*dNudy[0] + Ux[10]*dNudy[10] + Ux[11]*dNudy[11] + Ux[12]*dNudy[12] + Ux[13]*dNudy[13] + Ux[14]*dNudy[14] + Ux[15]*dNudy[15] + Ux[16]*dNudy[16] + Ux[17]*dNudy[17] + Ux[18]*dNudy[18] + Ux[19]*dNudy[19] + Ux[1]*dNudy[1] + Ux[20]*dNudy[20] + Ux[21]*dNudy[21] + Ux[22]*dNudy[22] + Ux[23]*dNudy[23] + Ux[24]*dNudy[24] + Ux[25]*dNudy[25] + Ux[26]*dNudy[26] + Ux[2]*dNudy[2] + Ux[3]*dNudy[3] + Ux[4]*dNudy[4] + Ux[5]*dNudy[5] + Ux[6]*dNudy[6] + Ux[7]*dNudy[7] + Ux[8]*dNudy[8] + Ux[9]*dNudy[9] + Uy[0]*dNudx[0] + Uy[10]*dNudx[10] + Uy[11]*dNudx[11] + Uy[12]*dNudx[12] + Uy[13]*dNudx[13] + Uy[14]*dNudx[14] + Uy[15]*dNudx[15] + Uy[16]*dNudx[16] + Uy[17]*dNudx[17] + Uy[18]*dNudx[18] + Uy[19]*dNudx[19] + Uy[1]*dNudx[1] + Uy[20]*dNudx[20] + Uy[21]*dNudx[21] + Uy[22]*dNudx[22] + Uy[23]*dNudx[23] + Uy[24]*dNudx[24] + Uy[25]*dNudx[25] + Uy[26]*dNudx[26] + Uy[2]*dNudx[2] + Uy[3]*dNudx[3] + Uy[4]*dNudx[4] + Uy[5]*dNudx[5] + Uy[6]*dNudx[6] + Uy[7]*dNudx[7] + Uy[8]*dNudx[8] + Uy[9]*dNudx[9];
  strain_rate_gp[4] = Ux[0]*dNudz[0] + Ux[10]*dNudz[10] + Ux[11]*dNudz[11] + Ux[12]*dNudz[12] + Ux[13]*dNudz[13] + Ux[14]*dNudz[14] + Ux[15]*dNudz[15] + Ux[16]*dNudz[16] + Ux[17]*dNudz[17] + Ux[18]*dNudz[18] + Ux[19]*dNudz[19] + Ux[1]*dNudz[1] + Ux[20]*dNudz[20] + Ux[21]*dNudz[21] + Ux[22]*dNudz[22] + Ux[23]*dNudz[23] + Ux[24]*dNudz[24] + Ux[25]*dNudz[25] + Ux[26]*dNudz[26] + Ux[2]*dNudz[2] + Ux[3]*dNudz[3] + Ux[4]*dNudz[4] + Ux[5]*dNudz[5] + Ux[6]*dNudz[6] + Ux[7]*dNudz[7] + Ux[8]*dNudz[8] + Ux[9]*dNudz[9] + Uz[0]*dNudx[0] + Uz[10]*dNudx[10] + Uz[11]*dNudx[11] + Uz[12]*dNudx[12] + Uz[13]*dNudx[13] + Uz[14]*dNudx[14] + Uz[15]*dNudx[15] + Uz[16]*dNudx[16] + Uz[17]*dNudx[17] + Uz[18]*dNudx[18] + Uz[19]*dNudx[19] + Uz[1]*dNudx[1] + Uz[20]*dNudx[20] + Uz[21]*dNudx[21] + Uz[22]*dNudx[22] + Uz[23]*dNudx[23] + Uz[24]*dNudx[24] + Uz[25]*dNudx[25] + Uz[26]*dNudx[26] + Uz[2]*dNudx[2] + Uz[3]*dNudx[3] + Uz[4]*dNudx[4] + Uz[5]*dNudx[5] + Uz[6]*dNudx[6] + Uz[7]*dNudx[7] + Uz[8]*dNudx[8] + Uz[9]*dNudx[9];
  strain_rate_gp[5] = Uy[0]*dNudz[0] + Uy[10]*dNudz[10] + Uy[11]*dNudz[11] + Uy[12]*dNudz[12] + Uy[13]*dNudz[13] + Uy[14]*dNudz[14] + Uy[15]*dNudz[15] + Uy[16]*dNudz[16] + Uy[17]*dNudz[17] + Uy[18]*dNudz[18] + Uy[19]*dNudz[19] + Uy[1]*dNudz[1] + Uy[20]*dNudz[20] + Uy[21]*dNudz[21] + Uy[22]*dNudz[22] + Uy[23]*dNudz[23] + Uy[24]*dNudz[24] + Uy[25]*dNudz[25] + Uy[26]*dNudz[26] + Uy[2]*dNudz[2] + Uy[3]*dNudz[3] + Uy[4]*dNudz[4] + Uy[5]*dNudz[5] + Uy[6]*dNudz[6] + Uy[7]*dNudz[7] + Uy[8]*dNudz[8] + Uy[9]*dNudz[9] + Uz[0]*dNudy[0] + Uz[10]*dNudy[10] + Uz[11]*dNudy[11] + Uz[12]*dNudy[12] + Uz[13]*dNudy[13] + Uz[14]*dNudy[14] + Uz[15]*dNudy[15] + Uz[16]*dNudy[16] + Uz[17]*dNudy[17] + Uz[18]*dNudy[18] + Uz[19]*dNudy[19] + Uz[1]*dNudy[1] + Uz[20]*dNudy[20] + Uz[21]*dNudy[21] + Uz[22]*dNudy[22] + Uz[23]*dNudy[23] + Uz[24]*dNudy[24] + Uz[25]*dNudy[25] + Uz[26]*dNudy[26] + Uz[2]*dNudy[2] + Uz[3]*dNudy[3] + Uz[4]*dNudy[4] + Uz[5]*dNudy[5] + Uz[6]*dNudy[6] + Uz[7]*dNudy[7] + Uz[8]*dNudy[8] + Uz[9]*dNudy[9];
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 2.0*eta_gp*strain_rate_gp[0];
  tau_gp[1] = 2.0*eta_gp*strain_rate_gp[1];
  tau_gp[2] = 2.0*eta_gp*strain_rate_gp[2];
  tau_gp[3] = eta_gp*strain_rate_gp[3];
  tau_gp[4] = eta_gp*strain_rate_gp[4];
  tau_gp[5] = eta_gp*strain_rate_gp[5];
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = strain_rate_gp[0] + strain_rate_gp[1] + strain_rate_gp[2];
  // y = A21.u at gauss point //
  Y[0] += -Np[0]*div_gp;
  Y[1] += -Np[1]*div_gp;
  Y[2] += -Np[2]*div_gp;
  Y[3] += -Np[3]*div_gp;
  
// total operations = 503
}
inline void MatMultMF_Stokes_MixedFEM3d_A22(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])
{
  const int nsd = 3;
  const int ntens = 6;
  double    p_gp;
  double    strain_rate_gp[6];
  double    tau_gp[6];
  double    div_gp;


  // pressure at gauss point //
  p_gp = Np[0]*P[0] + Np[1]*P[1] + Np[2]*P[2] + Np[3]*P[3];
  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //
  strain_rate_gp[0] = 0;
  strain_rate_gp[1] = 0;
  strain_rate_gp[2] = 0;
  strain_rate_gp[3] = 0;
  strain_rate_gp[4] = 0;
  strain_rate_gp[5] = 0;
  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //
  tau_gp[0] = 0;
  tau_gp[1] = 0;
  tau_gp[2] = 0;
  tau_gp[3] = 0;
  tau_gp[4] = 0;
  tau_gp[5] = 0;
  // y = A11.u at gauss point //
  // y = A12.p at gauss point //
  // divergence at gauss point //
  div_gp = 0;
  // y = A21.u at gauss point //
  Y[0] += 0;
  Y[1] += 0;
  Y[2] += 0;
  Y[3] += 0;
  
// total operations = 11
}
