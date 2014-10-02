
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
def test01():
	test_name = 'test01'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 2/ Check viscous block solve 
def test02():
	test_name = 'test02'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 3a/ Linear solve GalerkinMG
def test03a():
	test_name = 'test03a'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 3b/ Linear solve MFMG
def test03b():
	test_name = 'test03b'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 

# 3c/ Linear solve using MF,ASM,Galerkin and user specified level and direction dependent coarsening 
def test03c():
	test_name = 'test03c'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print('--------------------------------------------------------')
	print('-- Performing test:',test_name)
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 

# 3c-parallel/ Linear solve using MF,ASM,Galerkin and user specified level and direction dependent coarsening 
def test03cp(launcher):
	test_name = 'test03cp'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	ranks     = 8
	print('--------------------------------------------------------')
	print('-- Performing test:' , test_name , '[ processors = ', ranks ,']')

	space = ' '
	joblauncher = space.join( [launcher , '-np' , str(ranks)] )
	execcmd     = space.join( ['${PTATIN3D_DIR}/ptatin_driver_linear_ts.app' , '-options_file' ,  optsfile , '>' , ofile] )
	run         = space.join( [joblauncher , execcmd] )
	os.system(run)

	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print('---- Unix diff with expected file:',efile)
	os.system('diff ' + ofile+'.strip expected/expected.strip') 



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


	# set environment variable for bin directory
	if os.environ.get('PETSC_ARCH') == None:
		raise Exception('Environment variable PETSC_ARCH must be set to run unit tests')
	ptatin_bin_dir = os.path.join('..', '..', os.environ['PETSC_ARCH'], 'bin')
	os.environ['PTATIN3D_DIR'] = ptatin_bin_dir

	# determine path to mpiexec/mpirun
	launcher = determinempilaunchcommand()
	if launcher == None:
		has_mpi = False
	else:
		has_mpi = True


	# create output directories
	makedir('pt3dout')
	makedir('output')

	# remove contents of directory where all ptatin output is sent
	os.system('rm -f pt3dout/*')
	os.system('rm -f output/*')

	if test_name_target == 'all':
		test01()
		test02()
		test03a()
		test03b()
		test03c()
		if has_mpi == True:
			test03cp(launcher)

	elif test_name_target == 'test01':
		test01()

	elif test_name_target == 'test02':
		test02()

	elif test_name_target == 'test03a':
		test03a()

	elif test_name_target == 'test03b':
		test03b()

	elif test_name_target == 'test03c':
		test03c()

	elif test_name_target == 'test03cp':
		if has_mpi == True:
			test03cp(launcher)
	else:
		exit('-- Error: Unit test with name = "'+ test_name_target +'" was not found --')
	

if __name__ == "__main__":
	main()

