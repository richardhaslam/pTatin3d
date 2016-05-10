function nphase = Layer2tatin(prefix,filename,LayerPhase,LayerDepth)

nlayer = length(LayerPhase);

fid3 = fopen([filename,'.option'],'r+');
fseek(fid3, 0, 'eof');
fprintf(fid3,'#####################################\n');
fprintf(fid3,'#############LayerCake###############\n');
fprintf(fid3,'#####################################\n');

for i=0:nlayer-1
s1 = [prefix,'phaseLayer_',num2str(i)];
fprintf(fid3,'%s %d\n',s1,LayerPhase(i+1));
s1 = [prefix,'Ylayer_'];
fprintf(fid3,'%s %e\n',s1,LayerDepth(i+1));
end
fclose(fid3)
nphase = nlayer;
end