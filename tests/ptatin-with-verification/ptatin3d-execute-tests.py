
import sys
import os
import errno


def strip_keywords(infilename,outfilename):
	bad_words = ['(MB)', '(sec)']
	bad_words = bad_words + ['**']
	bad_words = bad_words + ['[pTatinModel]']
	bad_words = bad_words + ['[[']
	bad_words = bad_words + ['[pTatin]']

	#with open('t1.out') as oldfile, open('t1.strip.out', 'w') as newfile:

	oldfile = open(infilename,'r')
	newfile = open(outfilename,'w')

	for line in oldfile:
		if not any(bad_word in line for bad_word in bad_words):
			newfile.write(line)
			
	oldfile.close()
	newfile.close()


# Make a directory if it doesn't exist
def makedir(dirname):
	try:
		os.makedirs(dirname)
	except OSError as e:
		if e.errno == errno.EEXIST and os.path.isdir(dirname):
			pass
		else: raise

# Determine MPIEXEC
def determinempilaunchcommand():
	launch_cmd = None

	if os.environ.get('PETSC_DIR') == None:
		print('-- Warning: Environment variable PETSC_DIR is not set => all parallel tests will be ignored')
		return(launch_cmd)

	loc_petsc_var = os.path.join( os.environ['PETSC_DIR'] , os.environ['PETSC_ARCH'] , 'lib/petsc/conf/' , 'petscvariables' )
	#print(loc_petsc_var)

	file = open( loc_petsc_var, 'r' )
	for line in file:
		#Look for "MPIEXEC = " in petscvariables
		x = line.find('MPIEXEC =')
		if x != -1:
			broken_line = line.split()
			#print(broken_line)
			launch_cmd = broken_line[2]
			break
	file.close()
	if launch_cmd == '':
		exit('-- Error: Failed to determine how to launch parallel job --')
	
	return(launch_cmd)
	
#
# Test definitions
#

# 1/ Check mf operators are matching assembled operators
def test01(this_dir):
	test_name = 'test01'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	print('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')


# 2/ Check viscous block solve 
def test02(this_dir):
	test_name = 'test02'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	# Filter
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')


# 3a/ Linear solve GalerkinMG
def test03a(this_dir):
	test_name = 'test03a'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	# Filter
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')


# 3b/ Linear solve MFMG
def test03b(this_dir):
	test_name = 'test03b'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	# Filter
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')

# 3c/ Linear solve using MF,ASM,Galerkin and user specified level and direction dependent coarsening 
def test03c(this_dir):
	test_name = 'test03c'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ptatin_std_output_path + ' > ' + ofile)
	# Filter
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')

# 3c-parallel/ Linear solve using MF,ASM,Galerkin and user specified level and direction dependent coarsening 
def test03cp(this_dir,launcher):
	test_name = 'test03cp'
	ofile     = os.path.join(this_dir, 'output'   , test_name + '.output')
	efile     = os.path.join(this_dir, 'expected' , test_name + '.expected')
	optsfile  = os.path.join(this_dir, 'input'    , test_name + '.opts')
	ranks     = 8
	print('--------------------------------------------------------')
	print('-- Performing test:' , test_name , '[ processors = ', ranks ,']')

	space = ' '
	joblauncher = space.join( [launcher , '-np' , str(ranks)] )
	execcmd     = space.join( ['${PTATIN3D_DIR}/ptatin_driver_linear_ts.app' , '-options_file' ,  optsfile , ptatin_std_output_path, '>' , ofile] )
	run         = space.join( [joblauncher , execcmd] )
	os.system(run)

	# Filter
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,this_dir + '/expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip ' + this_dir + '/expected/expected.strip')



# 4/ SNES - MFMG



# Execute all tests
def main():
    
    
	if len(sys.argv) == 1:
		test_name_target = 'all'
		print('**')
		print('** Executing all pTatin3d unit tests')
		print('**')
	else:
		test_name_target = sys.argv[1]
		print('**')
		print('** Executing pTatin3d unit test: ' + test_name_target)
		print('**')


	if os.environ.get('PETSC_DIR') == None:
		raise Exception('Environment variable PETSC_DIR must be set to run unit tests')
	if os.environ.get('PETSC_ARCH') == None:
		raise Exception('Environment variable PETSC_ARCH must be set to run unit tests')

  # set environment variable for bin directory
#	ptatin_bin_dir = os.path.join('..', '..', os.environ['PETSC_ARCH'], 'bin')
#	os.environ['PTATIN3D_DIR'] = ptatin_bin_dir
#	print('ptatin_bin_dir',ptatin_bin_dir)

	path_to_this_file = os.path.realpath(__file__)
	#print(os.path.dirname(os.path.abspath(path_to_this_file)))
	this_dir = os.path.dirname(os.path.abspath(path_to_this_file))
	ptatin_bin_dir = os.path.join(this_dir,'..', '..', os.environ['PETSC_ARCH'], 'bin')
	os.environ['PTATIN3D_DIR'] = ptatin_bin_dir
	#print('[correct] ptatin_bin_dir',ptatin_bin_dir)

	os.environ['PTATIN3D_TEST_DIR'] = this_dir

	global ptatin_std_output_path
	ptatin_std_output_path = ' -output_path ' + os.path.join(this_dir,'pt3dout ')
	print(ptatin_std_output_path)

	# determine path to mpiexec/mpirun
	launcher = determinempilaunchcommand()
	if launcher == None:
		has_mpi = False
	else:
		has_mpi = True


	# create output directories
	o1_dir = os.path.join(this_dir,'pt3dout')
	o2_dir = os.path.join(this_dir,'output')
	makedir(o1_dir)
	makedir(o2_dir)

	# remove contents of directory where all ptatin output is sent
	os.system('rm -f ' + o1_dir + '/*')
	os.system('rm -f ' + o2_dir + '/*')

	if test_name_target == 'all':
		test01(this_dir)
		test02(this_dir)
		test03a(this_dir)
		test03b(this_dir)
		test03c(this_dir)
		if has_mpi == True:
			test03cp(this_dir,launcher)

	elif test_name_target == 'test01':
		test01(this_dir)

	elif test_name_target == 'test02':
		test02(this_dir)

	elif test_name_target == 'test03a':
		test03a(this_dir)

	elif test_name_target == 'test03b':
		test03b(this_dir)

	elif test_name_target == 'test03c':
		test03c(this_dir)

	elif test_name_target == 'test03cp':
		if has_mpi == True:
			test03cp(this_dir,launcher)
	else:
		exit('-- Error: Unit test with name = "'+ test_name_target +'" was not found --')
	

if __name__ == "__main__":
	main()

