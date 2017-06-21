echo "Here, look for all the differences to be small (<1e-13) and for errors"
echo 1:
srun -N 1 -n 1 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op cuda -mx 4 -my 4 -mz 4 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 2:
srun -N 1 -n 1 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op cuda -mx 13 -my 7 -mz 5 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 3:
srun -N 2 -n 2 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op cuda -mx 13 -my 17 -mz 5 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 4:
srun -N 1 -n 3 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 1 -my 1 -mz 3 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 5:
srun -N 1 -n 4 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 1 -my 2 -mz 8 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 6:
srun -N 1 -n 7 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 3 -my 7 -mz 11 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 7:
srun -N 2 -n 16 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 5 -my 7 -mz 14 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10 | grep A11_mfo

echo 8:
srun -N 2 -n 16 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 16 -my 16 -mz 8 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  -subrepart_frac 0.25 | grep A11_mfo

echo 9:
srun -N 2 -n 24 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 16 -my 16 -mz 16 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 10:
srun -N 6 -n 64 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 16 -my 16 -mz 16 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 11:
srun -N 6 -n 72 ../$PETSC_ARCH/bin/test_stokes_operators.app -compare_operators -a11_op subrepart -mx 128 -my 96 -mz 48 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10  | grep A11_mfo

echo 12:
srun -N 1 -n 1 ../$PETSC_ARCH/bin/test_stokes_operators.app -log_view -mx 16 -my 16 -mz 16 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10 -lattice_layout_perturb 0.0 -apply_A11mf_operator -iterations 2 -a11_op cuda > /dev/null

echo 13:
srun -N 1 -n 1 ../$PETSC_ARCH/bin/test_stokes_operators.app -log_view -mx 32 -my 32 -mz 32 -ptatin_model viscous_sinker -model_viscous_sinker_eta1 10 -lattice_layout_perturb 0.0 -apply_A11mf_operator -iterations 2 -a11_op subrepart > /dev/null

echo 14:
srun -N 1 -n 4 ../$PETSC_ARCH/bin/ptatin_driver_linear_ts.app -options_file ../src/models/viscous_sinker/examples/sinker-mfscaling.opts -a11_op avx -mx 16 -my 16 -mz 16 | grep "KSP Component"

echo 15:
srun -N 1 -n 1 ../$PETSC_ARCH/bin/ptatin_driver_linear_ts.app -options_file ../src/models/viscous_sinker/examples/sinker-mfscaling.opts -a11_op cuda -mx 16 -my 16 -mz 16 | grep "KSP Component"

echo 16:
srun -N 1 -n 12 ../$PETSC_ARCH/bin/ptatin_driver_linear_ts.app -options_file ../src/models/viscous_sinker/examples/sinker-mfscaling.opts -a11_op subrepart -mx 16 -my 16 -mz 16 | grep "KSP Component"
