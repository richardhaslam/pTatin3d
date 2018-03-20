

SUBROUTINE f_matmultmf_stokes_mixedfem3d_b11( FAC, eta_gp, Ux, Uy, Uz, P, Nu, dNudx, dNudy, dNudz, Np, Y )

  ! input
  DOUBLE PRECISION, INTENT(IN) :: FAC
  DOUBLE PRECISION, INTENT(IN) :: eta_gp

  DOUBLE PRECISION,  DIMENSION(0:80), INTENT(OUT) ::  Y
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  Ux
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  Uy
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  Uz
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  P
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  Nu
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  dNudx
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  dNudy
  DOUBLE PRECISION,  DIMENSION(0:26), INTENT(IN)  ::  dNudz
  DOUBLE PRECISION,  DIMENSION(0:3),  INTENT(IN)  ::  Np

  ! local variables
  DOUBLE PRECISION,  DIMENSION(0:5)  ::  strain_rate_gp
  DOUBLE PRECISION,  DIMENSION(0:5)  ::  tau_gp
  DOUBLE PRECISION  ::  div_gp


  ! strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point
  strain_rate_gp(0) = Ux(0)*dNudx(0) + Ux(10)*dNudx(10) + Ux(11)*dNudx(11) + Ux(12)*dNudx(12) + Ux(13)*dNudx(13) + Ux(14)*dNudx(14) + Ux(15)*dNudx(15) + Ux(16)*dNudx(16) + Ux(17)*dNudx(17) + Ux(18)*dNudx(18) + Ux(19)*dNudx(19) + Ux(1)*dNudx(1) + Ux(20)*dNudx(20) + Ux(21)*dNudx(21) + Ux(22)*dNudx(22) + Ux(23)*dNudx(23) + Ux(24)*dNudx(24) + Ux(25)*dNudx(25) + Ux(26)*dNudx(26) + Ux(2)*dNudx(2) + Ux(3)*dNudx(3) + Ux(4)*dNudx(4) + Ux(5)*dNudx(5) + Ux(6)*dNudx(6) + Ux(7)*dNudx(7) + Ux(8)*dNudx(8) + Ux(9)*dNudx(9)
  strain_rate_gp(1) = Uy(0)*dNudy(0) + Uy(10)*dNudy(10) + Uy(11)*dNudy(11) + Uy(12)*dNudy(12) + Uy(13)*dNudy(13) + Uy(14)*dNudy(14) + Uy(15)*dNudy(15) + Uy(16)*dNudy(16) + Uy(17)*dNudy(17) + Uy(18)*dNudy(18) + Uy(19)*dNudy(19) + Uy(1)*dNudy(1) + Uy(20)*dNudy(20) + Uy(21)*dNudy(21) + Uy(22)*dNudy(22) + Uy(23)*dNudy(23) + Uy(24)*dNudy(24) + Uy(25)*dNudy(25) + Uy(26)*dNudy(26) + Uy(2)*dNudy(2) + Uy(3)*dNudy(3) + Uy(4)*dNudy(4) + Uy(5)*dNudy(5) + Uy(6)*dNudy(6) + Uy(7)*dNudy(7) + Uy(8)*dNudy(8) + Uy(9)*dNudy(9)
  strain_rate_gp(2) = Uz(0)*dNudz(0) + Uz(10)*dNudz(10) + Uz(11)*dNudz(11) + Uz(12)*dNudz(12) + Uz(13)*dNudz(13) + Uz(14)*dNudz(14) + Uz(15)*dNudz(15) + Uz(16)*dNudz(16) + Uz(17)*dNudz(17) + Uz(18)*dNudz(18) + Uz(19)*dNudz(19) + Uz(1)*dNudz(1) + Uz(20)*dNudz(20) + Uz(21)*dNudz(21) + Uz(22)*dNudz(22) + Uz(23)*dNudz(23) + Uz(24)*dNudz(24) + Uz(25)*dNudz(25) + Uz(26)*dNudz(26) + Uz(2)*dNudz(2) + Uz(3)*dNudz(3) + Uz(4)*dNudz(4) + Uz(5)*dNudz(5) + Uz(6)*dNudz(6) + Uz(7)*dNudz(7) + Uz(8)*dNudz(8) + Uz(9)*dNudz(9)
  strain_rate_gp(3) = Ux(0)*dNudy(0) + Ux(10)*dNudy(10) + Ux(11)*dNudy(11) + Ux(12)*dNudy(12) + Ux(13)*dNudy(13) + Ux(14)*dNudy(14) + Ux(15)*dNudy(15) + Ux(16)*dNudy(16) + Ux(17)*dNudy(17) + Ux(18)*dNudy(18) + Ux(19)*dNudy(19) + Ux(1)*dNudy(1) + Ux(20)*dNudy(20) + Ux(21)*dNudy(21) + Ux(22)*dNudy(22) + Ux(23)*dNudy(23) + Ux(24)*dNudy(24) + Ux(25)*dNudy(25) + Ux(26)*dNudy(26) + Ux(2)*dNudy(2) + Ux(3)*dNudy(3) + Ux(4)*dNudy(4) + Ux(5)*dNudy(5) + Ux(6)*dNudy(6) + Ux(7)*dNudy(7) + Ux(8)*dNudy(8) + Ux(9)*dNudy(9) + Uy(0)*dNudx(0) + Uy(10)*dNudx(10) + Uy(11)*dNudx(11) + Uy(12)*dNudx(12) + Uy(13)*dNudx(13) + Uy(14)*dNudx(14) + Uy(15)*dNudx(15) + Uy(16)*dNudx(16) + Uy(17)*dNudx(17) + Uy(18)*dNudx(18) + Uy(19)*dNudx(19) + Uy(1)*dNudx(1) + Uy(20)*dNudx(20) + Uy(21)*dNudx(21) + Uy(22)*dNudx(22) + Uy(23)*dNudx(23) + Uy(24)*dNudx(24) + Uy(25)*dNudx(25) + Uy(26)*dNudx(26) + Uy(2)*dNudx(2) + Uy(3)*dNudx(3) + Uy(4)*dNudx(4) + Uy(5)*dNudx(5) + Uy(6)*dNudx(6) + Uy(7)*dNudx(7) + Uy(8)*dNudx(8) + Uy(9)*dNudx(9)
  strain_rate_gp(4) = Ux(0)*dNudz(0) + Ux(10)*dNudz(10) + Ux(11)*dNudz(11) + Ux(12)*dNudz(12) + Ux(13)*dNudz(13) + Ux(14)*dNudz(14) + Ux(15)*dNudz(15) + Ux(16)*dNudz(16) + Ux(17)*dNudz(17) + Ux(18)*dNudz(18) + Ux(19)*dNudz(19) + Ux(1)*dNudz(1) + Ux(20)*dNudz(20) + Ux(21)*dNudz(21) + Ux(22)*dNudz(22) + Ux(23)*dNudz(23) + Ux(24)*dNudz(24) + Ux(25)*dNudz(25) + Ux(26)*dNudz(26) + Ux(2)*dNudz(2) + Ux(3)*dNudz(3) + Ux(4)*dNudz(4) + Ux(5)*dNudz(5) + Ux(6)*dNudz(6) + Ux(7)*dNudz(7) + Ux(8)*dNudz(8) + Ux(9)*dNudz(9) + Uz(0)*dNudx(0) + Uz(10)*dNudx(10) + Uz(11)*dNudx(11) + Uz(12)*dNudx(12) + Uz(13)*dNudx(13) + Uz(14)*dNudx(14) + Uz(15)*dNudx(15) + Uz(16)*dNudx(16) + Uz(17)*dNudx(17) + Uz(18)*dNudx(18) + Uz(19)*dNudx(19) + Uz(1)*dNudx(1) + Uz(20)*dNudx(20) + Uz(21)*dNudx(21) + Uz(22)*dNudx(22) + Uz(23)*dNudx(23) + Uz(24)*dNudx(24) + Uz(25)*dNudx(25) + Uz(26)*dNudx(26) + Uz(2)*dNudx(2) + Uz(3)*dNudx(3) + Uz(4)*dNudx(4) + Uz(5)*dNudx(5) + Uz(6)*dNudx(6) + Uz(7)*dNudx(7) + Uz(8)*dNudx(8) + Uz(9)*dNudx(9)
  strain_rate_gp(5) = Uy(0)*dNudz(0) + Uy(10)*dNudz(10) + Uy(11)*dNudz(11) + Uy(12)*dNudz(12) + Uy(13)*dNudz(13) + Uy(14)*dNudz(14) + Uy(15)*dNudz(15) + Uy(16)*dNudz(16) + Uy(17)*dNudz(17) + Uy(18)*dNudz(18) + Uy(19)*dNudz(19) + Uy(1)*dNudz(1) + Uy(20)*dNudz(20) + Uy(21)*dNudz(21) + Uy(22)*dNudz(22) + Uy(23)*dNudz(23) + Uy(24)*dNudz(24) + Uy(25)*dNudz(25) + Uy(26)*dNudz(26) + Uy(2)*dNudz(2) + Uy(3)*dNudz(3) + Uy(4)*dNudz(4) + Uy(5)*dNudz(5) + Uy(6)*dNudz(6) + Uy(7)*dNudz(7) + Uy(8)*dNudz(8) + Uy(9)*dNudz(9) + Uz(0)*dNudy(0) + Uz(10)*dNudy(10) + Uz(11)*dNudy(11) + Uz(12)*dNudy(12) + Uz(13)*dNudy(13) + Uz(14)*dNudy(14) + Uz(15)*dNudy(15) + Uz(16)*dNudy(16) + Uz(17)*dNudy(17) + Uz(18)*dNudy(18) + Uz(19)*dNudy(19) + Uz(1)*dNudy(1) + Uz(20)*dNudy(20) + Uz(21)*dNudy(21) + Uz(22)*dNudy(22) + Uz(23)*dNudy(23) + Uz(24)*dNudy(24) + Uz(25)*dNudy(25) + Uz(26)*dNudy(26) + Uz(2)*dNudy(2) + Uz(3)*dNudy(3) + Uz(4)*dNudy(4) + Uz(5)*dNudy(5) + Uz(6)*dNudy(6) + Uz(7)*dNudy(7) + Uz(8)*dNudy(8) + Uz(9)*dNudy(9)

  ! deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point
  tau_gp(0) = 2.0*eta_gp*strain_rate_gp(0)
  tau_gp(1) = 2.0*eta_gp*strain_rate_gp(1)
  tau_gp(2) = 2.0*eta_gp*strain_rate_gp(2)
  tau_gp(3) = eta_gp*strain_rate_gp(3)
  tau_gp(4) = eta_gp*strain_rate_gp(4)
  tau_gp(5) = eta_gp*strain_rate_gp(5)

  ! y = A11.u at gauss point
  ! y = A12.p at gauss point
  ! divergence at gauss point
  div_gp = strain_rate_gp(0) + strain_rate_gp(1) + strain_rate_gp(2)

  ! y = A21.u at gauss point
  Y(0) = Y(0) + FAC*(dNudx(0)*tau_gp(0) + dNudy(0)*tau_gp(3) + dNudz(0)*tau_gp(4))
  Y(1) = Y(1) + FAC*(dNudx(0)*tau_gp(3) + dNudy(0)*tau_gp(1) + dNudz(0)*tau_gp(5))
  Y(2) = Y(2) + FAC*(dNudx(0)*tau_gp(4) + dNudy(0)*tau_gp(5) + dNudz(0)*tau_gp(2))
  Y(3) = Y(3) + FAC*(dNudx(1)*tau_gp(0) + dNudy(1)*tau_gp(3) + dNudz(1)*tau_gp(4))
  Y(4) = Y(4) + FAC*(dNudx(1)*tau_gp(3) + dNudy(1)*tau_gp(1) + dNudz(1)*tau_gp(5))
  Y(5) = Y(5) + FAC*(dNudx(1)*tau_gp(4) + dNudy(1)*tau_gp(5) + dNudz(1)*tau_gp(2))
  Y(6) = Y(6) + FAC*(dNudx(2)*tau_gp(0) + dNudy(2)*tau_gp(3) + dNudz(2)*tau_gp(4))
  Y(7) = Y(7) + FAC*(dNudx(2)*tau_gp(3) + dNudy(2)*tau_gp(1) + dNudz(2)*tau_gp(5))
  Y(8) = Y(8) + FAC*(dNudx(2)*tau_gp(4) + dNudy(2)*tau_gp(5) + dNudz(2)*tau_gp(2))
  Y(9) = Y(9) + FAC*(dNudx(3)*tau_gp(0) + dNudy(3)*tau_gp(3) + dNudz(3)*tau_gp(4))

  Y(10) = Y(10) + FAC*(dNudx(3)*tau_gp(3) + dNudy(3)*tau_gp(1) + dNudz(3)*tau_gp(5))
  Y(11) = Y(11) + FAC*(dNudx(3)*tau_gp(4) + dNudy(3)*tau_gp(5) + dNudz(3)*tau_gp(2))
  Y(12) = Y(12) + FAC*(dNudx(4)*tau_gp(0) + dNudy(4)*tau_gp(3) + dNudz(4)*tau_gp(4))
  Y(13) = Y(13) + FAC*(dNudx(4)*tau_gp(3) + dNudy(4)*tau_gp(1) + dNudz(4)*tau_gp(5))
  Y(14) = Y(14) + FAC*(dNudx(4)*tau_gp(4) + dNudy(4)*tau_gp(5) + dNudz(4)*tau_gp(2))
  Y(15) = Y(15) + FAC*(dNudx(5)*tau_gp(0) + dNudy(5)*tau_gp(3) + dNudz(5)*tau_gp(4))
  Y(16) = Y(16) + FAC*(dNudx(5)*tau_gp(3) + dNudy(5)*tau_gp(1) + dNudz(5)*tau_gp(5))
  Y(17) = Y(17) + FAC*(dNudx(5)*tau_gp(4) + dNudy(5)*tau_gp(5) + dNudz(5)*tau_gp(2))
  Y(18) = Y(18) + FAC*(dNudx(6)*tau_gp(0) + dNudy(6)*tau_gp(3) + dNudz(6)*tau_gp(4))
  Y(19) = Y(19) + FAC*(dNudx(6)*tau_gp(3) + dNudy(6)*tau_gp(1) + dNudz(6)*tau_gp(5))

  Y(20) = Y(20) + FAC*(dNudx(6)*tau_gp(4) + dNudy(6)*tau_gp(5) + dNudz(6)*tau_gp(2))
  Y(21) = Y(21) + FAC*(dNudx(7)*tau_gp(0) + dNudy(7)*tau_gp(3) + dNudz(7)*tau_gp(4))
  Y(22) = Y(22) + FAC*(dNudx(7)*tau_gp(3) + dNudy(7)*tau_gp(1) + dNudz(7)*tau_gp(5))
  Y(23) = Y(23) + FAC*(dNudx(7)*tau_gp(4) + dNudy(7)*tau_gp(5) + dNudz(7)*tau_gp(2))
  Y(24) = Y(24) + FAC*(dNudx(8)*tau_gp(0) + dNudy(8)*tau_gp(3) + dNudz(8)*tau_gp(4))
  Y(25) = Y(25) + FAC*(dNudx(8)*tau_gp(3) + dNudy(8)*tau_gp(1) + dNudz(8)*tau_gp(5))
  Y(26) = Y(26) + FAC*(dNudx(8)*tau_gp(4) + dNudy(8)*tau_gp(5) + dNudz(8)*tau_gp(2))
  Y(27) = Y(27) + FAC*(dNudx(9)*tau_gp(0) + dNudy(9)*tau_gp(3) + dNudz(9)*tau_gp(4))
  Y(28) = Y(28) + FAC*(dNudx(9)*tau_gp(3) + dNudy(9)*tau_gp(1) + dNudz(9)*tau_gp(5))
  Y(29) = Y(29) + FAC*(dNudx(9)*tau_gp(4) + dNudy(9)*tau_gp(5) + dNudz(9)*tau_gp(2))

  Y(30) = Y(30) + FAC*(dNudx(10)*tau_gp(0) + dNudy(10)*tau_gp(3) + dNudz(10)*tau_gp(4))
  Y(31) = Y(31) + FAC*(dNudx(10)*tau_gp(3) + dNudy(10)*tau_gp(1) + dNudz(10)*tau_gp(5))
  Y(32) = Y(32) + FAC*(dNudx(10)*tau_gp(4) + dNudy(10)*tau_gp(5) + dNudz(10)*tau_gp(2))
  Y(33) = Y(33) + FAC*(dNudx(11)*tau_gp(0) + dNudy(11)*tau_gp(3) + dNudz(11)*tau_gp(4))
  Y(34) = Y(34) + FAC*(dNudx(11)*tau_gp(3) + dNudy(11)*tau_gp(1) + dNudz(11)*tau_gp(5))
  Y(35) = Y(35) + FAC*(dNudx(11)*tau_gp(4) + dNudy(11)*tau_gp(5) + dNudz(11)*tau_gp(2))
  Y(36) = Y(36) + FAC*(dNudx(12)*tau_gp(0) + dNudy(12)*tau_gp(3) + dNudz(12)*tau_gp(4))
  Y(37) = Y(37) + FAC*(dNudx(12)*tau_gp(3) + dNudy(12)*tau_gp(1) + dNudz(12)*tau_gp(5))
  Y(38) = Y(38) + FAC*(dNudx(12)*tau_gp(4) + dNudy(12)*tau_gp(5) + dNudz(12)*tau_gp(2))
  Y(39) = Y(39) + FAC*(dNudx(13)*tau_gp(0) + dNudy(13)*tau_gp(3) + dNudz(13)*tau_gp(4))

  Y(40) = Y(40) + FAC*(dNudx(13)*tau_gp(3) + dNudy(13)*tau_gp(1) + dNudz(13)*tau_gp(5))
  Y(41) = Y(41) + FAC*(dNudx(13)*tau_gp(4) + dNudy(13)*tau_gp(5) + dNudz(13)*tau_gp(2))
  Y(42) = Y(42) + FAC*(dNudx(14)*tau_gp(0) + dNudy(14)*tau_gp(3) + dNudz(14)*tau_gp(4))
  Y(43) = Y(43) + FAC*(dNudx(14)*tau_gp(3) + dNudy(14)*tau_gp(1) + dNudz(14)*tau_gp(5))
  Y(44) = Y(44) + FAC*(dNudx(14)*tau_gp(4) + dNudy(14)*tau_gp(5) + dNudz(14)*tau_gp(2))
  Y(45) = Y(45) + FAC*(dNudx(15)*tau_gp(0) + dNudy(15)*tau_gp(3) + dNudz(15)*tau_gp(4))
  Y(46) = Y(46) + FAC*(dNudx(15)*tau_gp(3) + dNudy(15)*tau_gp(1) + dNudz(15)*tau_gp(5))
  Y(47) = Y(47) + FAC*(dNudx(15)*tau_gp(4) + dNudy(15)*tau_gp(5) + dNudz(15)*tau_gp(2))
  Y(48) = Y(48) + FAC*(dNudx(16)*tau_gp(0) + dNudy(16)*tau_gp(3) + dNudz(16)*tau_gp(4))
  Y(49) = Y(49) + FAC*(dNudx(16)*tau_gp(3) + dNudy(16)*tau_gp(1) + dNudz(16)*tau_gp(5))

  Y(50) = Y(50) + FAC*(dNudx(16)*tau_gp(4) + dNudy(16)*tau_gp(5) + dNudz(16)*tau_gp(2))
  Y(51) = Y(51) + FAC*(dNudx(17)*tau_gp(0) + dNudy(17)*tau_gp(3) + dNudz(17)*tau_gp(4))
  Y(52) = Y(52) + FAC*(dNudx(17)*tau_gp(3) + dNudy(17)*tau_gp(1) + dNudz(17)*tau_gp(5))
  Y(53) = Y(53) + FAC*(dNudx(17)*tau_gp(4) + dNudy(17)*tau_gp(5) + dNudz(17)*tau_gp(2))
  Y(54) = Y(54) + FAC*(dNudx(18)*tau_gp(0) + dNudy(18)*tau_gp(3) + dNudz(18)*tau_gp(4))
  Y(55) = Y(55) + FAC*(dNudx(18)*tau_gp(3) + dNudy(18)*tau_gp(1) + dNudz(18)*tau_gp(5))
  Y(56) = Y(56) + FAC*(dNudx(18)*tau_gp(4) + dNudy(18)*tau_gp(5) + dNudz(18)*tau_gp(2))
  Y(57) = Y(57) + FAC*(dNudx(19)*tau_gp(0) + dNudy(19)*tau_gp(3) + dNudz(19)*tau_gp(4))
  Y(58) = Y(58) + FAC*(dNudx(19)*tau_gp(3) + dNudy(19)*tau_gp(1) + dNudz(19)*tau_gp(5))
  Y(59) = Y(59) + FAC*(dNudx(19)*tau_gp(4) + dNudy(19)*tau_gp(5) + dNudz(19)*tau_gp(2))

  Y(60) = Y(60) + FAC*(dNudx(20)*tau_gp(0) + dNudy(20)*tau_gp(3) + dNudz(20)*tau_gp(4))
  Y(61) = Y(61) + FAC*(dNudx(20)*tau_gp(3) + dNudy(20)*tau_gp(1) + dNudz(20)*tau_gp(5))
  Y(62) = Y(62) + FAC*(dNudx(20)*tau_gp(4) + dNudy(20)*tau_gp(5) + dNudz(20)*tau_gp(2))
  Y(63) = Y(63) + FAC*(dNudx(21)*tau_gp(0) + dNudy(21)*tau_gp(3) + dNudz(21)*tau_gp(4))
  Y(64) = Y(64) + FAC*(dNudx(21)*tau_gp(3) + dNudy(21)*tau_gp(1) + dNudz(21)*tau_gp(5))
  Y(65) = Y(65) + FAC*(dNudx(21)*tau_gp(4) + dNudy(21)*tau_gp(5) + dNudz(21)*tau_gp(2))
  Y(66) = Y(66) + FAC*(dNudx(22)*tau_gp(0) + dNudy(22)*tau_gp(3) + dNudz(22)*tau_gp(4))
  Y(67) = Y(67) + FAC*(dNudx(22)*tau_gp(3) + dNudy(22)*tau_gp(1) + dNudz(22)*tau_gp(5))
  Y(68) = Y(68) + FAC*(dNudx(22)*tau_gp(4) + dNudy(22)*tau_gp(5) + dNudz(22)*tau_gp(2))
  Y(69) = Y(69) + FAC*(dNudx(23)*tau_gp(0) + dNudy(23)*tau_gp(3) + dNudz(23)*tau_gp(4))

  Y(70) = Y(70) + FAC*(dNudx(23)*tau_gp(3) + dNudy(23)*tau_gp(1) + dNudz(23)*tau_gp(5))
  Y(71) = Y(71) + FAC*(dNudx(23)*tau_gp(4) + dNudy(23)*tau_gp(5) + dNudz(23)*tau_gp(2))
  Y(72) = Y(72) + FAC*(dNudx(24)*tau_gp(0) + dNudy(24)*tau_gp(3) + dNudz(24)*tau_gp(4))
  Y(73) = Y(73) + FAC*(dNudx(24)*tau_gp(3) + dNudy(24)*tau_gp(1) + dNudz(24)*tau_gp(5))
  Y(74) = Y(74) + FAC*(dNudx(24)*tau_gp(4) + dNudy(24)*tau_gp(5) + dNudz(24)*tau_gp(2))
  Y(75) = Y(75) + FAC*(dNudx(25)*tau_gp(0) + dNudy(25)*tau_gp(3) + dNudz(25)*tau_gp(4))
  Y(76) = Y(76) + FAC*(dNudx(25)*tau_gp(3) + dNudy(25)*tau_gp(1) + dNudz(25)*tau_gp(5))
  Y(77) = Y(77) + FAC*(dNudx(25)*tau_gp(4) + dNudy(25)*tau_gp(5) + dNudz(25)*tau_gp(2))
  Y(78) = Y(78) + FAC*(dNudx(26)*tau_gp(0) + dNudy(26)*tau_gp(3) + dNudz(26)*tau_gp(4))
  Y(79) = Y(79) + FAC*(dNudx(26)*tau_gp(3) + dNudy(26)*tau_gp(1) + dNudz(26)*tau_gp(5))
  Y(80) = Y(80) + FAC*(dNudx(26)*tau_gp(4) + dNudy(26)*tau_gp(5) + dNudz(26)*tau_gp(2))

  ! total operations = 1058

END
