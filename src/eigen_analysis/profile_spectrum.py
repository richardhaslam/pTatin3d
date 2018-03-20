
import os
import math

def linspace(start, stop, n):
	if n == 1:
		yield stop
		return
	h = (stop - start) / (n - 1)
	for i in range(n):
		yield start + h * i


output_file = 'eigen_spectrum'
options_file = 'test_eigen_analysis_mfgmg.opts'
slepc_opts = '-eps_monitor_conv -eps_tol 1e-5'

job1_A = '-mx 16 -my 16 -mz 16 -dau_nlevels 3 -model_viscous_sinker_eta0 1.0 -fieldsplit_u_ksp_rtol 1.0e-1 -gbump_amp 0.0'
job2_A = '-mx 16 -my 16 -mz 16 -dau_nlevels 3 -model_viscous_sinker_eta0 1.0 -fieldsplit_u_ksp_rtol 1.0e-3 -gbump_amp 0.0'
job3_A = '-mx 16 -my 16 -mz 16 -dau_nlevels 3 -model_viscous_sinker_eta0 1.0 -fieldsplit_u_ksp_rtol 1.0e-6 -gbump_amp 0.0'

job1_B = '-mx 16 -my 16 -mz 16 -dau_nlevels 3 -model_viscous_sinker_eta0 1.0e-4 -fieldsplit_u_ksp_rtol 1.0e-1 -gbump_amp 0.0'
job2_B = '-mx 32 -my 32 -mz 32 -dau_nlevels 4 -model_viscous_sinker_eta0 1.0e-4 -fieldsplit_u_ksp_rtol 1.0e-3 -gbump_amp 0.0'
job3_B = '-mx 32 -my 32 -mz 32 -dau_nlevels 4 -model_viscous_sinker_eta0 1.0e-4 -fieldsplit_u_ksp_rtol 1.0e-6 -gbump_amp 0.0'



prefix = 'job1_A'  # INTPUT NAME
options_file = options_file + ' ' + job1_A # INPUT EXTRA OPTIONS

output_file = output_file + '_' + prefix


# Indicate which style of analysis required
analyser_opts = ' -eigen_anal 1'
# 1 - examine spectrum of stokes_operator.stokes_PC
# 2 - example spectrum of A11_operator.A11_PC
# 3 - example spectrum of A11smoother_operator.A11smoother_PC

# Set flag to indicate if the system is symmetric or not
is_symmetric = False

# Set flag to indicate if we want to look at Re/Im limits
sample_minmax  = False
sample_re_axis = False
sample_im_axis = False

sample_range   = True
nsample        = 4
sample_min_mag = 2.5e-1
sample_max_mag = 1.0


#  ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -n 2 ./ex1 -n 5000  -eps_monitor -eps_tol 1e-5  -eps_target_magnitude -eps_target 1.0e-9 -eps_nev 500 -eps_mpd 600

# ==========================================================================================
cmd_basic = './ptatin_eigen_analysis.app -options_file' + ' ' + options_file
cmd_basic = cmd_basic + ' ' + analyser_opts
cmd_basic = cmd_basic + ' ' + slepc_opts


# Examine boundaries of Re/Im axis
if sample_minmax == True:
	cmd = cmd_basic + ' -eps_nev 25 -eigen_anal_view_operators -eps_largest_magnitude > ' + output_file + '_mag_max'
	print 'Executing... ' + cmd
	ierr = os.system(cmd)

	cmd = cmd_basic + ' -eps_nev 25 -eps_smallest_magnitude > ' + output_file  + '_mag_min'
	print 'Executing... ' + cmd
	ierr = os.system(cmd)

if sample_re_axis == True:
	# Examine extents of Re axis
	cmd = cmd_basic + '-eps_nev 20 -eps_smallest_real > ' + output_file  + '_re_max'
	print 'Executing... ' + cmd
	ierr = os.system(cmd)

	cmd = cmd_basic + '-eps_nev 20 -eps_largest_real > ' + output_file + '_re_min'
	print 'Executing... ' + cmd
	ierr = os.system(cmd)


# Examine extents of Im axis
if sample_im_axis == True:
	if is_symmetric == False:
		cmd = cmd_basic + '-eps_nev 20 -eps_largest_imaginary > ' + output_file + '_im_max'
		print 'Executing... ' + cmd
		ierr = os.system(cmd)

		cmd = cmd_basic + '-eps_nev 20 -eps_smallest_imaginary > ' + output_file + '_im_min'
		print 'Executing... ' + cmd
		ierr = os.system(cmd)


# Examine range
if sample_range == True:
	l1 = math.log10(sample_min_mag)
	l2 = math.log10(sample_max_mag)

	list = linspace(l1,l2,nsample)
	for r in list:
		mag_to_target = math.pow(10.0,r)
		#print r,math.pow(10.0,r)

		cmd = cmd_basic + ' -eps_nev 10 -eps_target_magnitude -eps_target ' + str(mag_to_target) + ' > ' + output_file + '_mag_' + str(mag_to_target)
		print 'Executing... ' + cmd
		ierr = os.system(cmd)
