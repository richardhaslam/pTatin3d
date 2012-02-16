clear;
% parser for input file to be used with Gene3d model in pTatin3D
prefix = '-model_GENE3D_';
% How do you want to input the initial geometry? 0 => LayerCake 1=> extrude
% along Z from RGB file 2=> use a CAD file
initial_geom = 1; 
% what name of file do you want for the output ?
filename = 'gene3D_ext'
% what is the size of the model?
Origine = [0.0 0.0 0.0];
Length  = [12 5.5 2.0];
% what is resolution the model?
nel     = [24 12 4];
% what kind of rheological model are you using ? VISCOUS => 0 VISCOPLASTIC => 1
% VISCOPLASTIC WITH SOFTENING => 2
rheology_type = 1;
% for non linear models it is menditory to set viscosity cut off for the continuation
eta_min_cut_off = 1e-6;
eta_max_cut_off = 1.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the filename for the geometry if you have one
%fname = 'testcolor-01.tif';
fname = 'exemple.png';
%fname = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define your layering here if you want a layered cake
% enter phase index from bottom to top of the model
LayerPhase  = [1 2 1 2 1];
% enter the depth of the top interface of each layer
LayerDepth  = [1 2 3 4 6]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enter your rheological parameter for each lithology you wanna define
% there may be more phase index than lithologies each litho is a line
% Maybe this will be changed to a structure that looks more like a taras
% style library with lithologies being given a name instead of a line
% number I am thinking about it! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis.type    =    { '%e'    '%e'};
vis.options = { 'eta0_' 'rho_ '};
vis.val     =   [1      0 ; ...
                1e-3    0 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plas.type      = { '%e'    '%e'   '%e'   '%e' };
plas.options   = { 'C0_' 'Phi_ ' 'Tens_' 'Hs_'};
plas.val       = [ 0.27    0 0.01     1e+100 ; ...
                1e100   0 1e+100   1e+100 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soft.type    = { '%e'    '%e'           '%e'        '%e' };
soft.options = { 'Co_inf' 'Phi_inf_' 'eps_min_''eps_max_'};
soft.val     = [   0.027    0           0           0.1; ...
                   1e+100   0        1e+100         1e+100];
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF USER INTERFACE %%%%%%%%%%%%%%%%


%%%%%%%%% CREATE AND EXPORT THE OPT TEXT FILE %%%%%%%%%%%%%%%%%%%



fid3 = fopen([filename,'.option'],'w');
fprintf(fid3,'#####################################\n');
fprintf(fid3,'#############Generality###############\n');
fprintf(fid3,'#####################################\n');
s1 = [prefix,'rheol'];
fprintf(fid3,'%s %d\n',s1,rheology_type);

fprintf(fid3,'%s %d\n','-mx',nel(1));
fprintf(fid3,'%s %d\n','-my',nel(2));
fprintf(fid3,'%s %d\n','-mz',nel(3));

fprintf(fid3,'%s %d\n','-Lx',Length(1));
fprintf(fid3,'%s %d\n','-Ly',Length(2));
fprintf(fid3,'%s %d\n','-Lz',Length(3));

fprintf(fid3,'%s %d\n','-Ox',Origine(1));
fprintf(fid3,'%s %d\n','-Oy',Origine(2));
fprintf(fid3,'%s %d\n','-Oz',Origine(3));

if (rheology_type > 0 )
    s1 = [prefix,'eta_min_cut_off'];
    fprintf(fid3,'%s %e\n',s1,eta_min_cut_off);
    s1 = [prefix,'eta_max_cut_off'];
    fprintf(fid3,'%s %e\n',s1,eta_max_cut_off);
end

fclose(fid3);


switch (initial_geom)
    case 0
        nphase = Layer2tatin(prefix,filename,LayerPhase,LayerDepth);
    case 1
        LayerPhase=[];
        [nphase,LayerPhase] = RGB2tatin(filename,fname,LayerPhase,Origine(1),Origine(2),Length(1),Length(2));
    case 2
        error('not developped for now\n');
end

materialproperties_option_write(prefix,rheology_type,filename,nphase,vis,plas,soft,LayerPhase)

fid3 = fopen([filename,'.option'],'r+');
fseek(fid3, 0, 'eof');
fprintf(fid3,'########################################\n');
fprintf(fid3,'#User Happens numerics bellow this line#\n');
fprintf(fid3,'########################################\n');
fclose(fid3);

