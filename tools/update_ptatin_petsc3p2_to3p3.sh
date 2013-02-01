perl -p -i -e 's/private\/daimpl.h/petsc\-private\/daimpl.h/g' *.c
perl -p -i -e 's/private\/daimpl.h/petsc\-private\/daimpl.h/g' *.h

perl -p -i -e 's/private\/matimpl.h/petsc\-private\/matimpl.h/g' *.c
perl -p -i -e 's/private\/matimpl.h/petsc\-private\/matimpl.h/g' *.h


perl -p -i -e 's/DMGetInterpolation/DMCreateInterpolation/g' *.c
perl -p -i -e 's/DMGetMatrix/DMCreateMatrix/g' *.c
perl -p -i -e 's/DMGetInjection/DMCreateInjection/g' *.c
perl -p -i -e 's/DMGetInterpolationScale/DMCreateInterpolationScale/g' *.c
perl -p -i -e 's/PetscTypeCompare/PetscObjectTypeCompare/g' *.c

perl -p -i -e 's/ierr\s\=\sPCMGSetGalerkin\(pc\_i,PETSC\_FALSE\);CHKERRQ\(ierr\)\;/ierr = PCMGSetGalerkin(pc\_i,PETSC\_FALSE);CHKERRQ(ierr);\n\t\tierr = PCSetDM(pc_i,PETSC\_NULL);CHKERRQ(ierr);/g' *.c

perl -p -i -e 's/ierr\s=\sMatSetType\(\sA,\sMATAIJ\s\)\;CHKERRQ\(ierr\)\;/ierr = MatSetType( A, MATAIJ );CHKERRQ(ierr);\n\tierr = MatSetUp(A);CHKERRQ(ierr);/g' data_exchanger.c

perl -p -i -e 's/ierr\s=\sMatSetBlockSize\(B,3\)\;CHKERRQ\(ierr\)\;/\/*ierr = MatSetBlockSize(B,3);CHKERRQ(ierr);*\//g' stokes_operators.c

cd test_option_files
perl -p -i -e 's/da\_mat\_type/dm\_mat\_type/g' *.opts 
perl -p -i -e 's/schur\_factorization\_type/schur\_fact\_type/g' *.opts 
perl -p -i -e 's/chebychev/chebyshev/g' *.opts
