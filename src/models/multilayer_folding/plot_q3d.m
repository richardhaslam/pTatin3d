% plot
clear all, close all

% pseudo 2D cylindrical

lambda_x = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

q_th = [23.695022800 14.863101100 10.202873500 7.714419520 6.191435800 5.167653610 4.433314060 3.881244630 ...
        3.451220960 3.106864490];

q_comp = [23.869811300 14.731754700 10.071490800 7.611881460 6.154400640 5.154395960 4.411516220 3.864560300 ...
          3.500603100 3.059220390];

figure(1)
plot(lambda_x, q_th)
hold on
plot(lambda_x, q_comp, 'o')
xlabel('\lambda_x', 'fontsize', 14)
ylabel('q', 'fontsize', 14)
title('q vs \lambda_x for a pseudo 2D setup','fontsize', 18)

% compression in one direction and perturbation applied in both directions


H = 0.5e-2;

lambda_x = 0.05:0.01:1;
lambda_z = 0.05:0.01:1;

[lambda_x_2D, lambda_z_2D] = meshgrid(lambda_x, lambda_z);


kx = 2*pi./lambda_x_2D;
kz = 2*pi./lambda_z_2D;
lam = sqrt(kx.^2 + kz.^2);




mu_l = 100;                                                                   % viscosity of the layer
mu_m = 1;                                                                % viscosity of the matrix
exx  = -0.3;                                                                % strain rate in x
ezz  = -0.0;                                                                % strain rate in z
eyy  = -(exx+ezz);                                                          % strain rate in y assuming incompressibility

R = -exx/eyy;                                                               % normalized strain rate in x direction
v = mu_m/mu_l;                                                              % viscosity contrast
                                                                            % wavenumber
 kx_2 = kx.^2;
 kz_2 = kz.^2;
%
% k = (kx_2+kz_2).^0.5;
 k = lam*H;

 q = (-4*k)*(1-v)./ (2*k*(1-v.^2) -(1+v).^2*exp(k) + (1-v).^2*exp(-k));
 q3D_norm = - q/2.*( (kx_2./lam.^2)*exx + (kz_2./lam.^2)*ezz - eyy);



figure(2)
%[c,h] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm_2);


% theoritical growth rate
c = [7.01 4.43 3.71 3.07 2.71 2.37 2.24 2.04 1.92 1.70 1.52 1.24];

% computed growth rate
d = [6.90 4.34 3.74 3.01 2.70 2.35 2.24 2.05 1.93 1.74 1.51 1.23];

[c,h] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm,c,'LineColor','blue');
hold on
[d,j] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm,d,'--','LineColor','black');

xlabel('\lambda_x/H','fontsize', 14)
ylabel('\lambda_y/H','fontsize', 14)
title('growth rate q3d for R = 1','fontsize', 18)

clabel(c,h)


% compression in 2 directions

ezz2 = -0.3;
eyy2 = -(exx+ezz2);

q3D_norm_2 = - q/2.*( (kx_2./lam.^2)*exx + (kz_2./lam.^2)*ezz2 - eyy2);

figure(3)
%[c,h] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm_2);


% theoritical growth rate
e = [10.73 7.13 5.31 4.417 3.92 3.61 3.18];

% computed growth rate
f = [10.50 7.12 5.24 4.42 3.91 3.60 3.16];

[e,i] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm_2,e,'LineColor','blue');
hold on
[f,g] = contour(lambda_x_2D/H,lambda_z_2D/H,q3D_norm_2,f,'--','LineColor','black');

xlabel('\lambda_x/H','fontsize', 14)
ylabel('\lambda_y/H','fontsize', 14)
title('growth rate q3d for R = 1/2','fontsize', 18)

clabel(e,i)
