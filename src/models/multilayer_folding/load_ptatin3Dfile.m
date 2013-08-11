
function strucdata = load_ptatin3Dfile(dire,filename, m)
    addpath('/Users/marinec/LIB/petsc-3.2-p6/bin/matlab');
    data = PetscBinaryRead([dire, filename]);
    strucdata.x = data(1:3:length(data));
    strucdata.y = data(2:3:length(data));
    strucdata.z = data(3:3:length(data));

    strucdata.x = reshape(strucdata.x, m(1), m(2), m(3));
    strucdata.y = reshape(strucdata.y, m(1), m(2), m(3));
    strucdata.z = reshape(strucdata.z, m(1), m(2), m(3));
     