#!/bin/bash -l

#SBATCH --mail-user=dave.mayhem23@gmail.com
#SBATCH --mail-type=ALL

# Job Name (do not use a space!)
#SBATCH --job-name="pTat3d"

#SBATCH --output=p3d-%j.out
#SBATCH --error=p3d-%j.err

# number of mpi jobs (-n )
#SBATCH --ntasks=64
# number of cores per node - max=12 (-N ) 
#SBATCH --ntasks-per-node=8

# #SBATCH --cpu_bind=verbose,cores
#SBATCH --mem_bind=verbose,local

#SBATCH --time=00:30:00

#======START=====
module load slurm
echo "=================================================="
echo "  The current job ID is $SLURM_JOB_ID"
echo "  Running on $SLURM_NNODES nodes"
echo "  Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "  A total of $SLURM_NPROCS tasks is used"
echo "=================================================="


# no space
EXEC=${PWD}/ptatin_driver_asmsolve.app


aprun -n $SLURM_NPROCS -N $SLURM_NTASKS_PER_NODE $EXEC -ptatin_model viscous_sinker -mx 10 -my 10 -mz 10 -dau_nlevels 2 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type upper -stk_velocity_da_mat_type aij -ksp_type fgmres -fieldsplit_u_ksp_type fgmres -fieldsplit_u_mg_levels_ksp_max_it 2  -fieldsplit_u_ksp_monitor -ksp_monitor_true_residual -A11_operator_type 0,1,1,1,1 -model_viscous_sinker_eta0 1.0 -fieldsplit_u_mg_levels_pc_type jacobi -fieldsplit_u_mg_levels_ksp_max_it 1 -fieldsplit_p_ksp_type preonly -fieldsplit_p_pc_type jacobi -mx 96 -my 96 -mz 96  -fieldsplit_u_mg_coarse_pc_type lu -fieldsplit_u_mg_coarse_pc_factor_mat_solver_package superlu_dist    -da_processors_x 4 -da_processors_y 4 -da_processors_z 4 -dau_nlevels 5  -snes_view  -model_viscous_sinker_eta0 1.0e-6  -fieldsplit_u_mg_levels_ksp_type chebychev -fieldsplit_u_mg_levels_ksp_max_it 10 -fieldsplit_u_mg_levels_ksp_chebychev_estimate_eigenvalues 0,0.2,0,1.1 -ksp_max_it 10 -log_summary test.log


exit


