
import os

output_file = 'eigen_spectrum.dat'

options_file = 'test_basic_mfgmg.opts'

slepc_opts = '-eps_monitor_conv -eps_tol 1e-4 -eps_nev 20'

# Indicate which style of analysis required
analyser_opts = ' -eigen_anal 1'
# 1 - examine spectrum of stokes_operator.stokes_PC
# 2 - example spectrum of A11_operator.A11_PC
# 3 - example spectrum of A11smoother_operator.A11smoother_PC

# Set flag to indicate if the system is symmetric or not
is_symmetric = False

# Set flag to indicate if we want to look at Re/Im limits
examine_limits = False


#  ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec -n 2 ./ex1 -n 5000  -eps_monitor -eps_tol 1e-5  -eps_target_magnitude -eps_target 1.0e-9 -eps_nev 500 -eps_mpd 600

# ==========================================================================================
cmd_basic = './ptatin_eigen_analysis.app -options_file' + ' ' + options_file
cmd_basic = cmd_basic + ' ' + analyser_opts
cmd_basic = cmd_basic + ' ' + slepc_opts

print 'Basic launch command: ' + cmd_basic

# Examine boundaries of Re/Im axis
cmd = cmd_basic + ' -eps_nev 250 -eigen_anal_view_operators -eps_largest_magnitude > ' + output_file
print 'Executing... ' + cmd
ierr = os.system(cmd)

cmd = cmd_basic + ' -eps_nev 250 -eps_smallest_magnitude >> ' + output_file
print 'Executing... ' + cmd
#ierr = os.system(cmd)

if examine_limits == True:
	# Examine extents of Re axis
	cmd = cmd_basic + ' -eps_smallest_real >> ' + output_file
	print 'Executing... ' + cmd
	#ierr = os.system(cmd)

	cmd = cmd_basic + ' -eps_largest_real >> ' + output_file
	print 'Executing... ' + cmd
	#ierr = os.system(cmd)


	# Examine extents of Im axis
	if is_symmetric == False:
		cmd = cmd_basic + ' -eps_largest_imaginary >> ' + output_file
		print 'Executing... ' + cmd
		#ierr = os.system(cmd)

		cmd = cmd_basic + ' -eps_smallest_imaginary >> ' + output_file
		print 'Executing... ' + cmd
		#ierr = os.system(cmd)




