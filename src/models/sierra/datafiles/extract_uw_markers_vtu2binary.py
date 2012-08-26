
import sys
import random
import os.path
from optparse import OptionParser
from array import array

from paraview import servermanager
from paraview import vtk
from paraview import vtkConstants
from paraview.vtk import *
from paraview.vtk import io


# ===============================================================
def extract_coords_phase(ug_reader,outfilename):
	ug_markers = ug_reader.GetOutput()
	print '#  number of cells in mesh:', ug_markers.GetNumberOfCells()
	print '#  number of nodes = ', ug_markers.GetNumberOfPoints()
	if ug_markers.GetNumberOfCells() != ug_markers.GetNumberOfPoints():
		print '!!WARNING: This does not look like a marker file !!\n'
		#sys.exit(1)


	allpointdata = ug_markers.GetPointData()
	print allpointdata
	numpointdatatypes = allpointdata.GetNumberOfArrays()
	print numpointdatatypes
	for p in range(0,numpointdatatypes):
		print 'DataArrayName['+str(p)+']: ' + allpointdata.GetArrayName(p)
	
	

	pointdata = None
	pointdata = ug_markers.GetPointData().GetArray('phase')
	if pointdata == None:
		print '!!ERROR: Could not locate field with name "phase"'
		sys.exit(1)


	phase_list = []
	for c in range(0,ug_markers.GetNumberOfPoints()):
		phase_p = pointdata.GetValue(c)
		phase_list.append( int(phase_p) )


	coord_list = []
	nodes = ug_markers.GetPoints()
	for n in range(0,ug_markers.GetNumberOfPoints()):
		this_node = nodes.GetPoint(n)

		coord_list.append( this_node[0] )
		coord_list.append( this_node[1] )
		coord_list.append( this_node[2] )


	######################################
	name1 = 'coords_'+outfilename
	name2 = 'phase_'+outfilename
	
	print '# Writing file: ' + name1
	print '# Writing file: ' + name2

	binfile1 = open(name1,'wb')
	binfile2 = open(name2,'wb')

	# write header - number of points
	binfile1.write( str(ug_markers.GetNumberOfPoints()) + '\n' )
	binfile2.write( str(ug_markers.GetNumberOfPoints()) + '\n' )

	# write coords
	coord_array = array('d', coord_list)
	coord_array.tofile(binfile1)

	# write rock types
	phase_array = array('i', phase_list)
	phase_array.tofile(binfile2)

	binfile1.close()
	binfile2.close()



# ===============================================================
def extract_coords(ug_reader,outfilename):
	ug_markers = ug_reader.GetOutput()
	print '#  number of cells in mesh:', ug_markers.GetNumberOfCells()
	print '#  number of nodes = ', ug_markers.GetNumberOfPoints()
	if ug_markers.GetNumberOfCells() != ug_markers.GetNumberOfPoints():
		print '!!WARNING: This does not look like a marker file !!\n'
		#sys.exit(1)

	coord_list = []
	nodes = ug_markers.GetPoints()
	for n in range(0,ug_markers.GetNumberOfPoints()):
		this_node = nodes.GetPoint(n)

		coord_list.append( this_node[0] )
		coord_list.append( this_node[1] )
		coord_list.append( this_node[2] )


	######################################
	name = 'coords_' + outfilename
	
	print '# Writing file: ' + name

	binfile = open(name,'wb')

	# write vtk data type index: FORCE AS DOUBLE
	binfile.write( '11 \n' )

	# write header - number of points
	binfile.write( str(ug_markers.GetNumberOfPoints()) + '\n' )

	# write coords
	coord_array = array('d', coord_list)
	coord_array.tofile(binfile)

	binfile.close()


# DataType index
# bit, unsigned_char, char, unsigned_short, short, unsigned_int, int, unsigned_long, long, float, double
##
## TAKEN FROM vtkType.h
# #define VTK_VOID 0
# #define VTK_BIT 1
# #define VTK_CHAR 2
# #define VTK_SIGNED_CHAR 15
# #define VTK_UNSIGNED_CHAR 3
# #define VTK_SHORT 4
# #define VTK_UNSIGNED_SHORT 5
# #define VTK_INT 6
# #define VTK_UNSIGNED_INT 7
# #define VTK_LONG 8
# #define VTK_UNSIGNED_LONG 9
# #define VTK_FLOAT 10
# #define VTK_DOUBLE 11
# #define VTK_ID_TYPE 12
##
##
# 6 => int (Int32)
# 8 => long int (Int64)
#
# 11 => double (Float64)
#
# ===============================================================
def extract_field(ug_reader,fieldname,outfilename):
	ug_markers = ug_reader.GetOutput()
	if ug_markers.GetNumberOfCells() != ug_markers.GetNumberOfPoints():
		print '!!WARNING: This does not look like a marker file !!\n'
		#sys.exit(1)

	allpointdata = ug_markers.GetPointData()
	numpointdatatypes = allpointdata.GetNumberOfArrays()
	print '#  Number of data arrays: '+ str(numpointdatatypes)
	for p in range(0,numpointdatatypes):
		print '#    DataArrayName['+str(p)+']: ' + allpointdata.GetArrayName(p)
	
	print '# Looking for datafield name: ' + fieldname
	

	pointdata = None
	pointdata = ug_markers.GetPointData().GetArray(fieldname)
#	print pointdata
	if pointdata == None:
		print '!!ERROR: Could not locate field with name "' + str(fieldname) + '"'
		sys.exit(1)
	else:
		print '# Found datafield name: ' + fieldname
		print '# Data type index: ' + str( pointdata.GetDataType() )

	

	field_list = []
	for c in range(0,ug_markers.GetNumberOfPoints()):
		field_p = pointdata.GetValue(c)
#		field_list.append( int(field_p) )
#		print field_p
		field_list.append( field_p )


	######################################
	name = fieldname + '_' + outfilename
	
	print '# Writing file: ' + name

	binfile = open(name,'wb')

	datatypeindex = pointdata.GetDataType()

# #define VTK_VOID 0
# #define VTK_BIT 1
# #define VTK_CHAR 2
# #define VTK_SIGNED_CHAR 15
# #define VTK_UNSIGNED_CHAR 3
# #define VTK_SHORT 4
# #define VTK_UNSIGNED_SHORT 5
# #define VTK_INT 6
# #define VTK_UNSIGNED_INT 7
# #define VTK_LONG 8
# #define VTK_UNSIGNED_LONG 9
# #define VTK_FLOAT 10
# #define VTK_DOUBLE 11
# #define VTK_ID_TYPE 12

	# write field types
	if datatypeindex   == 0:  # VTK_VOID
		print '!!ERROR: Unsupported type VTK_VOID'
		sys.exit(1)

	elif datatypeindex == 1:  # VTK_BIT
		print '!!ERROR: Unsupported type VTK_BIT'
		sys.exit(1)

	elif datatypeindex == 2:  # VTK_CHAR
		field_array = array('c', field_list)	

	elif datatypeindex == 3:  # VTK_SIGNED_CHAR
		field_array = array('b', field_list)	

	elif datatypeindex == 4:  # VTK_SHORT
		field_array = array('h', field_list)	

	elif datatypeindex == 5:  # VTK_UNSIGNED_SHORT
		field_array = array('H', field_list)	

	elif datatypeindex == 6:  # VTK_INT
		field_array = array('i', field_list)	

	elif datatypeindex == 7:  # VTK_UNSIGNED_INT
		field_array = array('I', field_list)	

	elif datatypeindex == 8:  # VTK_LONG
		field_array = array('l', field_list)	

	elif datatypeindex == 9:  # VTK_UNSIGNED_LONG
		field_array = array('L', field_list)	

	elif datatypeindex == 10: # VTK_FLOAT
		field_array = array('f', field_list)	

	elif datatypeindex == 11: # VTK_DOUBLE
		field_array = array('d', field_list)

	elif datatypeindex == 12: # VTK_ID_TYPE
		print '!!ERROR: Unsupported type VTK_ID_TYPE'
		sys.exit(1)

	else:
		print '!!ERROR: Unsupported data type: ' + str(datatypeindex)
		sys.exit(1)

	# write vtk data type index
	binfile.write( str(datatypeindex) + '\n' )
	
	# write header - number of points
	binfile.write( str(ug_markers.GetNumberOfPoints()) + '\n' )

	# write data to file	
	field_array.tofile(binfile)

	binfile.close()


# ===============================================================
def extact_field(ug_markers,fieldname,outfilename):

	pointdata = None
	pointdata = ug_markers.GetPointData().GetArray(fieldname)
#	print pointdata
	if pointdata == None:
		print '!!ERROR: Could not locate field with name "' + str(fieldname) + '"'
		sys.exit(1)
	else:
		print '  # Found datafield name: ' + fieldname
		print '  # Data type index: ' + str( pointdata.GetDataType() )

	field_list = []
	for c in range(0,ug_markers.GetNumberOfPoints()):
		field_p = pointdata.GetValue(c)
		field_list.append( field_p )


	######################################
	name = fieldname + '_' + outfilename
	
	print '  # Writing file: ' + name

	binfile = open(name,'wb')

	datatypeindex = pointdata.GetDataType()

	# write field types
	if datatypeindex   == 0:  # VTK_VOID
		print '!!ERROR: Unsupported type VTK_VOID'
		sys.exit(1)

	elif datatypeindex == 1:  # VTK_BIT
		print '!!ERROR: Unsupported type VTK_BIT'
		sys.exit(1)

	elif datatypeindex == 2:  # VTK_CHAR
		field_array = array('c', field_list)	

	elif datatypeindex == 3:  # VTK_SIGNED_CHAR
		field_array = array('b', field_list)	

	elif datatypeindex == 4:  # VTK_SHORT
		field_array = array('h', field_list)	

	elif datatypeindex == 5:  # VTK_UNSIGNED_SHORT
		field_array = array('H', field_list)	

	elif datatypeindex == 6:  # VTK_INT
		field_array = array('i', field_list)	

	elif datatypeindex == 7:  # VTK_UNSIGNED_INT
		field_array = array('I', field_list)	

	elif datatypeindex == 8:  # VTK_LONG
		field_array = array('l', field_list)	

	elif datatypeindex == 9:  # VTK_UNSIGNED_LONG
		field_array = array('L', field_list)	

	elif datatypeindex == 10: # VTK_FLOAT
		field_array = array('f', field_list)	

	elif datatypeindex == 11: # VTK_DOUBLE
		field_array = array('d', field_list)

	elif datatypeindex == 12: # VTK_ID_TYPE
		print '!!ERROR: Unsupported type VTK_ID_TYPE'
		sys.exit(1)

	else:
		print '!!ERROR: Unsupported data type: ' + str(datatypeindex)
		sys.exit(1)

	# write vtk data type index
	binfile.write( str(datatypeindex) + '\n' )
	
	# write header - number of points
	binfile.write( str(ug_markers.GetNumberOfPoints()) + '\n' )

	# write data to file	
	field_array.tofile(binfile)

	binfile.close()



# ===============================================================
def extract_field_prompt(ug_reader,outfilename):

	ug_markers = ug_reader.GetOutput()
	if ug_markers.GetNumberOfCells() != ug_markers.GetNumberOfPoints():
		print '!!WARNING: This does not look like a marker file !!\n'
		#sys.exit(1)

	allpointdata = ug_markers.GetPointData()
	numpointdatatypes = allpointdata.GetNumberOfArrays()
	print '#  Number of data arrays: '+ str(numpointdatatypes)
	for p in range(0,numpointdatatypes):
		print '#    DataArrayName['+str(p)+']: ' + allpointdata.GetArrayName(p)
	print ''
		
	fieldname = None
	while True:
		fieldname = raw_input("  Select a field to extract... (or exit to quit) ")
		if fieldname == 'exit':
			print '  Finished extraction'
			break
		else:
			extact_field(ug_markers,fieldname,outfilename)





def main():

	# create a built-in connection
	if not servermanager.ActiveConnection:
			connection = servermanager.Connect()


	optparser=OptionParser(usage='usage: %prog -i <filename1> -o <filename2>',
												 add_help_option=True,
												 description="""Read particle vtu file """ + 
												 """and write data to a single binary file.""")

	optparser.add_option( "-i", "--input", dest="opt_inputfile",
										help="Input file name", metavar="FILE")

	optparser.add_option( "-o", "--output", dest="opt_outputfile",
										help="Output file name", metavar="FILE")

	optparser.add_option( "-f", "--fields", dest="opt_fieldlist",
										help="Field list to extract", metavar="FILE")


	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
			optparser.print_help()
			sys.exit(1)

	infilename = options.opt_inputfile
	print '# Reading file: ' + infilename

	reader = io.vtkXMLUnstructuredGridReader()
	reader.SetFileName(infilename)
	reader.Update()

	outfilename = options.opt_outputfile
	if os.path.splitext(infilename)[1]=='.vtu':
		outfilename = os.path.splitext(infilename)[0] + "_bin.dat"

	
# DEPRECIATED
#	extract_coords_phase(reader,options.opt_outputfile)

	# DUMP COORDS
	extract_coords(reader,outfilename)

	# DUMP FIELDS VIA COMMAND LINE ARGS
#	extract_field(reader,'phase',outfilename)
#	extract_field(reader,'wil',outfilename)
#	extract_field(reader,'eta',outfilename)
	if options.opt_fieldlist != None:
		#print options.opt_fieldlist
		#flist = os.path.splitext(options.opt_fieldlist)==','
		flist = options.opt_fieldlist.split(',')
		#print flist
		for f in flist:
			print '========================================================================'
			extract_field(reader,str(f),outfilename)

	else:
		# DUMP FIELDS VIA PROMPT
		print '========================================================================'
		extract_field_prompt(reader,outfilename)

		
if __name__ == '__main__':
	main()


