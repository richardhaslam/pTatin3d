function materialproperties_option_write(prefix,rheology_type,filename,nphase,vis,plas,soft,i_litho)

switch rheology_type
    case 0
        val     = vis.val ;
        options = vis.options;
        type    = vis.type;
    case 1
        val     = [vis.val,plas.val];
        options = [vis.options,plas.options];
        type    = [vis.type,plas.type];
    case 2
        val     = [vis.val,plas.val,soft.val];
        options = [vis.options,plas.options,soft.options];
        type    = [vis.type,plas.type,soft.type];
end

fid3 = fopen([filename,'.option'],'r+');
fseek(fid3, 0, 'eof');
fprintf(fid3,'#####################################\n');
fprintf(fid3,'#############Materials###############\n');
fprintf(fid3,'#####################################\n');

s1 = [prefix,'nphase'];
fprintf(fid3,'%s %d\n',s1,nphase);

for i_phase = 0:nphase-1
    fprintf(fid3,'#### param for phase indexed %d #####\n',i_phase);
    for i_opt = 1:size(val,2)
        s = [prefix,char(options(i_opt)),num2str(i_phase)];
        fprintf(fid3,['%s ',char(type(i_opt)),'\n'],s,val(i_litho(i_phase+1),i_opt));
    end
    fprintf(fid3,'#####################################\n');
end
fclose(fid3);
end