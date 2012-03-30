#!/bin/sh

# @ job_name = pTat3s
# @ comment = "pTat3d_profile"

# ## Blue Gene job
# @ job_type = bluegene

# @ error = $(job_name)_$(jobid).err
# @ output = $(job_name)_$(jobid).out
# @ environment = COPY_ALL;

# ## Syntax: wall_clock_limit = <hardlimit,softlimit>
# ##   [[hours:]minutes:]seconds[.fraction]
# @ wall_clock_limit = 0:20:00

# ## Syntax: notification = always|error|start|never|complete
# ##   where:
# ##   - always : Notify the user when the job begins, ends,
# ##     or if it incurs error conditions.
# ##   - error : Notify the user only if the job fails.
# ##   - start : Notify the user only when the job begins.
# ##   - never : Never notify the user.
# ##   - complete : Notify the user only when the job ends.
# ## Default value: complete
# @ notification = always
# @ notify_user = dave.mayhem23@gmail.com

# ## Syntax:  class = short|prod|mid
# @ class = short

# ## Syntax: bg_size = number_of_nodes
# @ bg_size = 1024
# ## Syntax: bg_connection = TORUS | PREFER_TORUS | MESH
# @ bg_connection = TORUS

# @ queue
mpirun -cwd /bgscratch/$USER/dmay -mode VN ${PWD}/ptatin_driver_asmsolve.app -options_file ${PWD}/optsfile -ptatin_model viscous_sinker -mx 10 -my 10 -mz 10 -dau_nlevels 2 -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_factorization_type upper -stk_velocity_da_mat_type aij -ksp_type fgmres -fieldsplit_u_ksp_type fgmres -fieldsplit_u_mg_levels_ksp_max_it 2  -fieldsplit_u_ksp_monitor -ksp_monitor_true_residual -A11_operator_type 0,1,1,1,1 -model_viscous_sinker_eta0 1.0 -fieldsplit_u_mg_levels_pc_type jacobi -fieldsplit_u_mg_levels_ksp_max_it 1 -fieldsplit_p_ksp_type preonly -fieldsplit_p_pc_type jacobi -mx 96 -my 96 -mz 96  -fieldsplit_u_mg_coarse_pc_type lu -fieldsplit_u_mg_coarse_pc_factor_mat_solver_package superlu_dist    -da_processors_x 4 -da_processors_y 4 -da_processors_z 4 -dau_nlevels 5  -snes_view  -model_viscous_sinker_eta0 1.0e-6  -fieldsplit_u_mg_levels_ksp_type chebychev -fieldsplit_u_mg_levels_ksp_max_it 10 -fieldsplit_u_mg_levels_ksp_chebychev_estimate_eigenvalues 0,0.2,0,1.1 -ksp_max_it 10 -log_summary test.log


