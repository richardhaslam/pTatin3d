#pushd src

git ls-files '*.[ch]' | xargs perl -pi -e '
  s,(\(\(PetscObject\)[^)]+\))->comm,PetscObjectComm$1,;
  s,\(\(PetscObject\)\(([^()]+)\)\)->comm,PetscObjectComm((PetscObject)$1),;
  s,\bPETSC_NULL\b,NULL,g;
  s,DMDABoundaryType,DMBoundaryType,g;
  s,DMDA_BOUNDARY_,DM_BOUNDARY_,g;
  s,DMDAGetCoordinateDA,DMGetCoordinateDM,;
  s,DMDAGetGhostedCoordinates,DMGetCoordinatesLocal,;
  s,DMDAGetCoordinates,DMGetCoordinates,;
  s,petsc-private/daimpl.h,petsc-private/dmdaimpl.h,;
  s/(DMDASetUniformCoordinates.*) NULL,NULL/$1 0.,0./;
  s/(DMDASetUniformCoordinates.*) NULL,NULL/$1 0.,0./;
  s/PetscSynchronizedFlush\( comm \)/PetscSynchronizedFlush( comm, PETSC_STDOUT )/;
  s,([^-])private/(snes|ksp|pc),$1petsc-private/$2,;
  s,snes->xtol,snes->stol,;
  s,SNES_CONVERGED_PNORM_RELATIVE,SNES_CONVERGED_SNORM_RELATIVE,;
  s@PetscMalloc\(([^,;]*[^,; ]) *\* *sizeof\([^,;()]+\),@PetscMalloc1($1,@;
  s@PetscMalloc2\(([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+)\)@PetscMalloc2($1,$3,$4,$6)@;
  s@PetscMalloc3\(([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+),([^,;()]+)\)@PetscMalloc3($1,$3,$4,$6,$7,$9)@;
  s@PetscNew\([^,;()]+ *, *@PetscNew(@;
  s@PetscNewLog\(([^,;()]+) *,[^,;()]+, *@PetscNewLog($1,@;
  s@KSPDefaultGetWork@KSPSetWorkVecs@;
  s@LogObjectMemory\((ksp)@LogObjectMemory((PetscObject)$1@;
  s@LogObjectParent\((pc)@LogObjectParent((PetscObject)$1@;
  s@PetscObjectComposeFunctionDynamic\(([^,]+),([^,]+),[^,]+,([^;]+)\);@PetscObjectComposeFunction($1,$2,$3);@;
  s@KSPSkipConverged@KSPConvergedSkip@;
  s@VecCreateMPIWithArray\(([^,]+)@VecCreateMPIWithArray($1,1@;
  s@PetscObject(Take|Grant)Access@PetscObjectSAWs$1Access@;
  s@KSPDefaultBuildSolution@KSPBuildSolutionDefault@;
  s@KSPDefaultBuildResidual@KSPBuildResidualDefault@;
  s@KSPDefaultDestroy@KSPDestroyDefault@;
  s@KSPRegisterDynamic\(([^,]+),[^,]+,[^,]+,([^,]+)\)@KSPRegister($1,$2)@;
  s@PCRegisterDynamic\(([^,]+),[^,]+,[^,]+,([^,]+)\)@PCRegister($1,$2)@;
  s@^(\s+)ierr = DMCreateMatrix\(([^,]+),([^,]+),@$1ierr = DMSetMatType($2,$3);CHKERRQ(ierr);\n$1ierr = DMCreateMatrix($2,@;
  s@PetscGetTime\(@PetscTime\(@;
  s@MatGetArray@MatSeqAIJGetArray@;
  s@MatRestoreArray@MatSeqAIJRestoreArray@;
  s@^.*DMView_DA_Private.*\n@@;
  s@((?:KSP|PC)[SG]etOperators)\(([^,;]+),([^,;]+),([^,;]+),([^,;]+)\)@$1($2,$3,$4)@;
'

git apply - <<"EOF"
diff --git i/src/dmda_redundant.c w/src/dmda_redundant.c
index 6956580..04c48a0 100644
--- i/src/dmda_redundant.c
+++ w/src/dmda_redundant.c
@@ -183,12 +183,11 @@ PetscErrorCode DMDACreate3dRedundant(DM da,PetscInt si, PetscInt ei, PetscInt sj
 	ierr = x_DMDACreate3d( PETSC_COMM_SELF, wrapNP, st, (ei-si),(ej-sj),(ek-sk), PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE, n_dofs, sw, 0,0,0, &_sda );CHKERRQ(ierr);
 	/* more hacky shit due to DMSetFromOptions() in constructor */
 	{
-		DM_DA          *dd = (DM_DA*)_sda->data;
 		DM dd_da_coordinates;
 		PetscInt _m,_n,_p,_M,_N,_P;
 		ierr = DMDAGetInfo(_sda,0,&_m,&_n,&_p,&_M,&_N,&_P,0,0,0,0,0,0);CHKERRQ(ierr);
 		ierr = x_DMDACreate3d(PetscObjectComm((PetscObject)_sda),wrapNP,st,_m,_n,_p,_M,_N,_P,3,sw,0,0,0,&dd_da_coordinates);CHKERRQ(ierr);
-		dd->da_coordinates = dd_da_coordinates;
+		ierr = DMSetCoordinateDM(_sda,dd_da_coordinates);CHKERRQ(ierr);
 	
 	}
 	ierr = DMDASetUniformCoordinates( _sda, 0.0,1.0, 0.0,1.0, 0.0,1.0 );CHKERRQ(ierr);
EOF
