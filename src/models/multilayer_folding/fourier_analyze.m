clear all
close all
addpath('/Users/marinec/LIB/petsc-3.2-p6/bin/matlab');


model = 'test_multilayer_ptatin3d';
type = 'lx50lz50';

scale_length = 100;
H = 0.5;
el_res_i = 16;
el_res_j = [6,2,6];
el_res_k = 16;

% nx = ;
% ny = ;
% nz = ;

el_resolution = [sum(el_res_i),sum(el_res_j),sum(el_res_k)];
node_resolution = 2*el_resolution + 1;


coord_filename = 'ptat3dcpf.dmda-velocity-coords_step000001';
vel_filename = 'ptat3dcpf.Xu_step000001';
dire = ['/Users/marinec/',model,'/',type,'/'];
%dire = ['/Users/marinec/ptatin3d/src/test1/']

%[mesh.x, mesh.y, mesh.z] =
%mesh = struct('x',[],'y', [], 'z',[]);
mesh = load_ptatin3Dfile(dire, coord_filename, node_resolution);

mesh.x = scale_length*mesh.x;
mesh.y = scale_length*mesh.y;
mesh.z = scale_length*mesh.z;


%slice of the bottom layer
slice = 2;
sj = [ 2*el_res_j(1) + 1, ...
       2*(el_res_j(1) + el_res_j(2)) + 1 ];

figure(1)

N = node_resolution(1);
M = node_resolution(3);

X_ = reshape(mesh.x(:,sj(1),:), N,M);
Y_ = reshape(mesh.y(:,sj(1),:), N,M);
Z_ = reshape(mesh.z(:,sj(1),:), N,M);

figure(1)
surf(X_,Z_,Y_)
axis equal;
xlabel('x')
ylabel('Z')


YY    = Y_-50;
A0    = H/100;
time  = 5.0000e-05;
time_ = time*scale_length;
eyy   = 0.4;

A = (max(max(YY))-min(min(YY)))/2
q3D  = (log(A/A0)/time_) / eyy
figure(2)
surf(X_,Z_,YY)

%AA = A0 * exp (26 * eyy* time )


% xmax = max(max(mesh.x(:,sj(1),1)));
% xmin = min(min(mesh.x(:,sj(1),1)));
% zmax = max(max(mesh.z(1,sj(1),:)));
% zmin = min(min(mesh.z(1,sj(1),:)));
%
% dx = (xmax-xmin)/(N-1);
% dz = (zmax-zmin)/(M-1);
%
% xs = xmin:dx:xmax;
% new_mesh.x = repmat(xs',1,M);
%
% zs = zmin:dz:zmax;
% new_mesh.z = repmat(zs,N,1);


%meshz = mesh.z(:,sj(1),:) - (max(max(mesh.z(:,sj(1),:))) + min(min(mesh.z(:,sj(1),:))))/2;

%--------------------Need a regular domain sampling for the FFT--------------

% F = TriScatteredInterp(reshape(mesh.x(:,sj(1),:),M*N,1),reshape(mesh.z(:,sj(1),:),M*N,1),reshape(meshz,M*N,1));
% new_mesh.y = F(new_mesh.x,new_mesh.z);
%
% %remove border: Could be fucked up due to the interpolation (NaN)
%
% new_mesh.x = new_mesh.x(2:N-1,2:M-1);
% new_mesh.y = new_mesh.y(2:N-1,2:M-1);
% new_mesh.z = new_mesh.z(2:N-1,2:M-1);
% figure(2)
% surf(new_mesh.x, new_mesh.z, new_mesh.y);
% axis equal;
%
% %--------------------FFT-----------------------
%
%
%
%
% zmax = max(max(new_mesh.z));
% zmin = min(min(new_mesh.z));
%
% nN = N-2;
% nM = M-2;
%
% NFFT = 2^nextpow2(nN);
%
% ndz = (zmax-zmin)/(nM-1);
%
% meshz_fft = repmat(zmin:ndz:zmax,NFFT,1);
%
% fe = 1/dx;
% f=fe*(0:NFFT-1)/NFFT;
% l = (0:NFFT-1)/NFFT;
%
% freq = repmat(f',1,nM);
% lambda = repmat(l',1,nM);
%
%
% ffty= abs(fft(new_mesh.y,NFFT,1));
% figure(3)
% surf(1./freq, meshz_fft, ffty)
