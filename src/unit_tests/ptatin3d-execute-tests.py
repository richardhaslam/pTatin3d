
import sys
import os


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


#
# Test definitions
#

# 1/ Check mf operators are matching assembled operators
def test01():
	test_name = 'test01'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print '--------------------------------------------------------'
	print '-- Performing test:',test_name
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print '---- Unix diff with expected file:',efile
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 2/ Check viscous block solve 
def test02():
	test_name = 'test02'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print '--------------------------------------------------------'
	print '-- Performing test:',test_name
	os.system('${PTATIN3D_DIR}/test_stokes_operators.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print '---- Unix diff with expected file:',efile
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 3a/ Linear solve GalerkinMG
def test03a():
	test_name = 'test03a'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print '--------------------------------------------------------'
	print '-- Performing test:',test_name
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print '---- Unix diff with expected file:',efile
	os.system('diff ' + ofile+'.strip expected/expected.strip') 


# 3b/ Linear solve MFMG
def test03b():
	test_name = 'test03b'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print '--------------------------------------------------------'
	print '-- Performing test:',test_name
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print '---- Unix diff with expected file:',efile
	os.system('diff ' + ofile+'.strip expected/expected.strip') 

# 3c/ Linear solve using MF,ASM,Galerkin and user specified level and direction dependent coarsening 
def test03c():
	test_name = 'test03c'
	ofile     = './output/'   + test_name + '.output'
	efile     = './expected/' + test_name + '.expected'
	optsfile  = './input/'    + test_name + '.opts'
	print '--------------------------------------------------------'
	print '-- Performing test:',test_name
	os.system('${PTATIN3D_DIR}/ptatin_driver_linear_ts.app -options_file ' + optsfile + ' > ' + ofile)
	# Filter 
	strip_keywords(ofile,ofile+'.strip')
	strip_keywords(efile,'expected/expected.strip')
	# DIFFERENCES
	print '---- Unix diff with expected file:',efile
	os.system('diff ' + ofile+'.strip expected/expected.strip') 



# 4/ SNES - MFMG



# Execute all tests
def main():
	print '**'
	print '** Executing pTatin3d unit tests'
	print '**'

	# set environment variable for bin directory
	ptatin_bin_dir = os.path.join('..', '..', os.environ['PETSC_ARCH'], 'bin')
	os.environ['PTATIN3D_DIR'] = ptatin_bin_dir

	# remove contents of directory where all ptatin output is sent
	os.system('rm -f pt3dout/*')
	os.system('rm -f output/*')

	test01()
	test02()
	test03a()
	test03b()
	test03c()


if __name__ == "__main__":
	main()

