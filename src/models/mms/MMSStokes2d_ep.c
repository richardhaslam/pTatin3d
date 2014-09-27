#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_u"
PetscErrorCode MMSStokes2D_u(PetscReal x,PetscReal y,PetscReal u[])
{
  u[0] = pow(x, 3)*y + pow(x, 2) + x*y + x;
  u[1] = -1.5*pow(x, 2)*pow(y, 2) - 2.0*x*y - 0.5*pow(y, 2) - y;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_p"
PetscErrorCode MMSStokes2D_p(PetscReal x,PetscReal y,PetscReal p[])
{
  p[0] = pow(x, 2)*pow(y, 2) + x*y - 0.361111111111111;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_gradu"
PetscErrorCode MMSStokes2D_gradu(PetscReal x,PetscReal y,PetscReal gradu[])
{
  gradu[0] = 3*pow(x, 2)*y + 2*x + y + 1;
  gradu[1] = pow(x, 3) + x;
  gradu[2] = -3.0*x*pow(y, 2) - 2.0*y;
  gradu[3] = -3.0*pow(x, 2)*y - 2.0*x - 1.0*y - 1;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_eta"
PetscErrorCode MMSStokes2D_eta(PetscReal x,PetscReal y,PetscReal eta[])
{
  eta[0] = -1.0*sin(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + 1.001;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_Fu"
PetscErrorCode MMSStokes2D_Fu(PetscReal x,PetscReal y,PetscReal fu[])
{
  fu[0] = -2*x*pow(y, 2) - y + (-3.0*x*y - 1.0)*(-2.0*sin(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + 2.002) + (6*x*y + 2)*(-2.0*sin(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + 2.002) - 2.0*(2.0*x*pow(y, 2) + 1.0*y)*(3*pow(x, 2)*y + 2*x + y + 1)*cos(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) - 2.0*(2.0*pow(x, 2)*y + 1.0*x)*(0.5*pow(x, 3) - 1.5*x*pow(y, 2) + 0.5*x - 1.0*y)*cos(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111);
  fu[1] = -2*pow(x, 2)*y - x + (-3.0*pow(x, 2) - 1.0)*(-2.0*sin(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + 2.002) - 2.0*(2.0*x*pow(y, 2) + 1.0*y)*(0.5*pow(x, 3) - 1.5*x*pow(y, 2) + 0.5*x - 1.0*y)*cos(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) - 2.0*(2.0*pow(x, 2)*y + 1.0*x)*(-3.0*pow(x, 2)*y - 2.0*x - 1.0*y - 1)*cos(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + (-2.0*sin(1.0*pow(x, 2)*pow(y, 2) + 1.0*x*y - 0.361111111111111) + 2.002)*(1.5*pow(x, 2) - 1.5*pow(y, 2) + 0.5);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MMSStokes2D_Fp"
PetscErrorCode MMSStokes2D_Fp(PetscReal x,PetscReal y,PetscReal fp[])
{
  fp[0] = 0;
  PetscFunctionReturn(0);
}

