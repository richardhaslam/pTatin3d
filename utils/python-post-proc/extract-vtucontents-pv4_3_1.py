
# http://stackoverflow.com/questions/17854468/paraview-python-scripting-equivalent-of-file-save-data
# http://www.cfd-online.com/Forums/paraview/96308-extracting-paraview-data-into-python-arrays.html


import os
import sys
from optparse import OptionParser

from paraview.numpy_support import vtk_to_numpy
from paraview import simple


def ExtractVTUContents(inputfile):

	tracer_ug = simple.XMLUnstructuredGridReader( FileName=inputfile )
	tracer_ug.UpdatePipeline()
	tracer_ug.UpdatePipelineInformation()


	rawdata = simple.servermanager.Fetch(tracer_ug)


	import vtk.numpy_interface.dataset_adapter as dsa

	# Wrap the raw data object to access NumPy friendly API
	rawdata = dsa.WrapDataObject(rawdata)

	dataInfo = tracer_ug.GetDataInformation()
	print('[PV file type] ' + dataInfo.GetDataSetTypeAsString() + '\n')

	pointDataInfo = dataInfo.GetPointDataInformation()
	print('[PV file contents]')
	print(pointDataInfo)
	print('\n')


	pa = pointDataInfo.GetArrayInformation(0)          # access by index
	pa = pointDataInfo.GetArrayInformation("pressure") # OR access by name
	#print(pa)

	pointData = tracer_ug.PointData
	#print("[pointData methods]")
	#print(dir(pointData))
	#print('\n')

	print("[Point Data Fields]")
	for n in range(pointData.GetNumberOfArrays()):
		print('  Name: ' + pointData.GetArray(n).GetName() + '; Data Range ' +
		str(pointData.GetArray(n).GetRange()))
	print('\n')


	#print( data.GetPointData().GetArray("pressure").GetValue(1) )
	field = rawdata.GetPointData().GetArray("pressure")
	print('Number of pressure vals: ' + str(rawdata.GetPointData().GetArray("pressure").GetSize()))


	#print("[Field methods]")
	#print dir(field)


	pressurefield = rawdata.PointData["pressure"]

	if isinstance(pressurefield, dsa.VTKNoneArray):
		print("Warning: Field being extracted does not exist");
	else:
		print('Pressure value at index 10: ' + str(pressurefield[10]))



def main():

	# parse command line options
	optparser=OptionParser(usage='usage: %prog -i <filename1>',
			add_help_option=True,
			description="""Example of reading a VTU file """ +
			"""and extracting one field into a numpy array.""")

	optparser.add_option( "-i", "--input", dest="opt_inputfile",
			help="Input file name", metavar="FILE")

	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
			optparser.print_help()
			sys.exit(1)

	infilename = options.opt_inputfile

	if os.path.splitext(infilename)[1]=='.vtu':
		ExtractVTUContents(infilename)
	else:
		print('Warning: Input file specified is not a valid VTU file')
		print('Warning: A valid VTU file must have the extension .vtu')
		print('Warning: Found input file: ' + infilename)
		optparser.print_help()
		sys.exit(1)


if __name__ == "__main__":
	main()
