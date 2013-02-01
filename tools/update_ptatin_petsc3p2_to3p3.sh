perl -p -i -e 's/private\/daimpl.h/petsc\-private\/daimpl.h/g' *.c
perl -p -i -e 's/private\/daimpl.h/petsc\-private\/daimpl.h/g' *.h

perl -p -i -e 's/private\/matimpl.h/petsc\-private\/matimpl.h/g' *.c
perl -p -i -e 's/private\/matimpl.h/petsc\-private\/matimpl.h/g' *.h


perl -p -i -e 's/DMGetInterpolation/DMCreateInterpolation/g' *.c
perl -p -i -e 's/DMGetMatrix/DMCreateMatrix/g' *.c
perl -p -i -e 's/DMGetInjection/DMCreateInjection/g' *.c
perl -p -i -e 's/DMGetInterpolationScale/DMCreateInterpolationScale/g' *.c
perl -p -i -e 's/PetscTypeCompare/PetscObjectTypeCompare/g' *.c

cd test_option_files
perl -p -i -e 's/schur\_factorization\_type/schur\_fact\_type/g' *.opts 

