
import sys
from optparse import OptionParser


from paraview.servermanager import *


def convert_pvtu2vtu(infilename,outfilename):
	# Unstructured Grid (particles)
	reader   = sources.XMLPUnstructuredGridReader(FileName=infilename)

	# DataMode = {0,1,2} => {ascii, encoded, appended}
	writer   = writers.XMLUnstructuredGridWriter(Input=reader,DataMode=1,FileName=outfilename)

	writer.UpdatePipeline()


def main():
	Connect()

	optparser=OptionParser(usage='usage: %prog -i <filename1>',
												 add_help_option=True,
												 description="""Read particle vtu file """ + 
												 """and write data to a single binary file.""")

	optparser.add_option( "-i", "--input", dest="opt_inputfile",
										help="Input file name", metavar="FILE")

	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
			optparser.print_help()
			sys.exit(1)


	infilename = options.opt_inputfile
	
	if os.path.splitext(infilename)[1]=='.pvtu':
		vts_name = os.path.splitext(infilename)[0] + ".vtu"

		print 'Reading: ' + infilename
		print 'Writing: ' + vts_name
		convert_pvtu2vtu(infilename,vts_name)


if __name__ == '__main__':
	main()
