from os import environ
import ptatin3d_execute_tests

# Generate batch files 
# Note that these produce output on a different (faster) filesystem
# But currently the main output still goes into a subfolder here on the home system. TODO is fixing that.
def generate_batch():
  testNames = ['test01','test02','test03a','test03b','test03c','test03cp']
  outPutDir = '/scratch/daint/' + environ['USER']
  PETSC_ARCH = environ['PETSC_ARCH'] #error if not defined
  for testName in testNames :
    outPutFile = outPutDir + '/'+testName+'.out'
    f = open(testName + '.batch','w')
    f.write('#!/bin/bash\n')
    if testName == 'test03cp':
      nprocstr='8'
    else :
      nprocstr = '1'
    f.write('#SBATCH --ntasks='+nprocstr+'\n')
    f.write('#SBATCH --time=00:05:00\n')
    f.write('#SBATCH --output='+outPutFile +'\n')
    if testName == 'test01' or testName == 'test02' :
      appName = 'test_stokes_operators.app'
    else :
      appName = 'ptatin_driver_linear_ts.app'
    f.write('aprun -n ' + nprocstr +' ../../'+PETSC_ARCH+'/bin/' + appName + ' -options_file ./input/' + testName + '.opts\n')
    f.close()

    # Generate a script to run and check output
    g = open(testName + '_check.py','w')
    g.write('#A script to compare the output of  ' + testName + '.batch and compare the output to the reference\n');
    g.write('import os\n')
    g.write('from ptatin3d_execute_tests import strip_keywords\n')
    g.write('def check_' + testName + '() :\n')
    eFile = './expected/' + testName + '.expected'
    g.write('  strip_keywords(\"' + outPutFile + '\",\"' + outPutFile + '.strip\")\n')
    g.write('  strip_keywords(\"' + eFile      + '\",\"' + eFile      + '.strip\")\n')
    g.write('  print(\'---- Unix diff with expected file ' + eFile + ' ---\')\n')
    g.write('  os.system(\'diff ' + outPutFile +'.strip ./expected/' + testName +'.expected.strip\')\n')
    g.write('if __name__==\"__main__\":\n')
    g.write('  check_' + testName + '()\n')
    g.close()
  print('To perform one of the tests, \n 1. sbatch testXXX.batch\n 2. python XXX_check.py\n')
  print('where XXX is an element of\n')
  print str(testNames) + '\n'

if __name__=="__main__":
  generate_batch()
