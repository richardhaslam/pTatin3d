
import os
import sys
from optparse import OptionParser


from paraview import simple
#from paraview.servermanager import *


def convert_pvtu2vtu(infilename,outfilename):
	# Unstructured Grid (particles)
	reader   = simple.XMLPartitionedUnstructuredGridReader(FileName=infilename)

	# DataMode = {0,1,2} => {ascii, encoded, appended}
	writer   = simple.XMLUnstructuredGridWriter(Input=reader,DataMode=1,FileName=outfilename)

	writer.UpdatePipeline()


def main():
	#Connect()

	optparser=OptionParser(usage='usage: %prog -i <filename1>',
			add_help_option=True,
			description="""Read parallel vtu (PVTU) file """ +
			"""and write data to a single binary VTU file.""")

	optparser.add_option( "-i", "--input", dest="opt_inputfile",
			help="Input file name", metavar="FILE")

	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
			optparser.print_help()
			sys.exit(1)


	infilename = options.opt_inputfile

	if os.path.splitext(infilename)[1]=='.pvtu':
		vts_name = os.path.splitext(infilename)[0] + ".vtu"

		print('Reading: ' + infilename)
		print('Writing: ' + vts_name)
		convert_pvtu2vtu(infilename,vts_name)
	else:
		print('Warning: Input file specified is not a valid PVTU file')
		print('Warning: A valid PVTU file must have the extension .pvtu')
		print('Warning: Found input file: ' + infilename)
		optparser.print_help()
		sys.exit(1)


if __name__ == '__main__':
  pvV = simple.GetParaViewVersion()
  majorV = pvV.major
  minorV = pvV.minor
  if majorV < 4:
    print('VersionError: Script is not tested for ParaView versions < v4.3')
    sys.exit(1)
  elif majorV == 4:
    if minorV < 3:
      print('VersionError: Script is not tested for ParaView versions < v4.3')
      sys.exit(1)

  main()
