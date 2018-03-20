function [nphase,i_litho] = RGB2tatin(filename,fname,i_litho,x0,y0,xend,yend)
% the file must be  RGB file, desactivate the anti-aliasing while saving
% it works fine with *.tif and *.png
% do not use jpg unless you want a thousands of phases
% the rest you can try
% on Illustrator save the file using artbox and create one that fits the
% model size
if isempty(fname)
    error('User must provide a valid filename for the image')
end
I=imread(fname);
figure(1);image(I);
pause(1);
[X,map] = rgb2ind(I(1:end-1,1:end-1,:),100);
nphase = size(map,1);
fprintf(1, 'number of phase found %d \n', nphase);
[nx,ny]=size(X);
for i=1:nx
    X2(nx-i+1,:)=X(i,:);
end
clear('X','I');
figure(2);
pcolor(double(X2));shading flat; colorbar; title('phase number for the input file');
X2=X2';
[nx,ny]=size(X2);
if isempty(i_litho)
    i_litho = input(' Enter the Array \n that contains the lithology index \n for each phase index as they appear ordered on the picture\n');
end

%%%%%%%%% EXPORT THE PMAP FILE %%%%%%%%%%%%%%%%%%%
fid2 = fopen([filename,'.pmap'],'w');
fprintf(fid2,'number vary first in x and than y \n');
fprintf(fid2,'%d\n', nx);
fprintf(fid2,'%d\n' ,ny);
fprintf(fid2,'%1.8e %1.8e %1.8e %1.8e\n',x0,y0,xend,yend);
fprintf(fid2,'%d ',X2);
fclose(fid2);