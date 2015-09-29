
# http://stackoverflow.com/questions/17854468/paraview-python-scripting-equivalent-of-file-save-data
# http://www.cfd-online.com/Forums/paraview/96308-extracting-paraview-data-into-python-arrays.html


import os
import sys
from optparse import OptionParser
import numpy as np
from paraview.numpy_support import vtk_to_numpy
from paraview import simple

import vtk.numpy_interface.dataset_adapter as dsa


def ExtractDefGridAsVTU(inputfile,inputbasefile,outputfile):
	
	tracer_ug = simple.XMLUnstructuredGridReader( FileName=inputfile )
	tracer_ug.UpdatePipeline()
	tracer_ug.UpdatePipelineInformation()

	rawdata_tracer = simple.servermanager.Fetch(tracer_ug)

	# Wrap the raw data object to access NumPy friendly API
	rawdata_tracer = dsa.WrapDataObject(rawdata_tracer)


	indexfield = rawdata_tracer.PointData["index"]
	#indexfield = indexfield -1
	npoints_tracer = rawdata_tracer.GetNumberOfPoints()

	defmesh_ug = simple.XMLStructuredGridReader( FileName= inputbasefile )
	defmesh_ug.UpdatePipeline()
	defmesh_ug.UpdatePipelineInformation()

	rawdata_defmesh = simple.servermanager.Fetch(defmesh_ug)
	rawdata_defmesh = dsa.WrapDataObject(rawdata_defmesh)

	npoints_defmesh = rawdata_defmesh.GetNumberOfPoints() 

	print('# tracers from input vtu file: ' + str(npoints_tracer))
	print('# mesh nodes from base vts file: ' + str(npoints_defmesh))

	if npoints_tracer != npoints_defmesh:
		print('ERROR: Particle dataset does not contain the same number of points as the base VTS file')
		sys.exit(1)

	coord_t = rawdata_tracer.GetPoints()
	coord_dm = rawdata_defmesh.GetPoints()

	for ii in range(0, npoints_tracer):
		tracer_idx = indexfield[ii]
		coord_dm[tracer_idx,0] = coord_t[ii,0]
		coord_dm[tracer_idx,1] = coord_t[ii,1]
		coord_dm[tracer_idx,2] = coord_t[ii,2]

	writer = simple.XMLStructuredGridWriter(Input=defmesh_ug,DataMode=1,FileName=outputfile)
	writer.UpdatePipeline()



def main():
	
	# parse command line options
	optparser=OptionParser(usage='usage: %prog -i <filename1> -b <filename2> -o <filename3>',
			add_help_option=True,
			description="""Example of reading a VTU file """ + 
			"""and extracting a deformed mesh.""")

	optparser.add_option( 	"-i", "--input", dest="opt_inputfile",
				help="Input vtu file name", metavar="FILE")

	optparser.add_option( 	"-b", "--inputbase", dest="opt_inputbasefile",
				help="Input base vts file name", metavar="FILE")

	optparser.add_option( 	"-o", "--output", dest="opt_outputfile",
				help="Output vts file name", metavar="FILE")

	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
		optparser.print_help()
		sys.exit(1)

	if options.opt_inputbasefile == None:
		optparser.print_help()
		sys.exit(1)

	if options.opt_outputfile == None:
		optparser.print_help()
		sys.exit(1)
	
	infilename     = options.opt_inputfile
	inbasefilename = options.opt_inputbasefile
	outputfilename = options.opt_outputfile

	if os.path.splitext(infilename)[1] != '.vtu':
		print('Warning: Input file specified is not a valid VTU file')
		print('Warning: A valid VTU file must have the extension .vtu')
		print('Warning: Found input file: ' + infilename)
		optparser.print_help()
		sys.exit(1)

	if os.path.splitext(inbasefilename)[1] != '.vts':
		print('Warning: Input base file specified is not a valid VTS file')
		print('Warning: A valid VTS file must have the extension .vts')
		print('Warning: Found input file: ' + inbasefilename)
		optparser.print_help()
		sys.exit(1)

	if os.path.splitext(outputfilename)[1] != '.vts':
		print('Warning: Output file specified is not valid')
		print('Warning: A valid VTS file must have the extension .vts')
		print('Warning: Found input file: ' + outputfilename)
		optparser.print_help()
		sys.exit(1)
	

	print('[Reading VTU file: ' + infilename + ']')
	print('[Reading VTS base file: ' + inbasefilename + ']')
	print('[Writing VTS file: ' + outputfilename + ']')
	ExtractDefGridAsVTU(infilename,inbasefilename,outputfilename)


if __name__ == "__main__":
	main()
