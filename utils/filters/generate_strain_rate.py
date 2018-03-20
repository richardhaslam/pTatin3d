
import sys
#import numpy
from optparse import OptionParser

from paraview.servermanager import *

from paraview import vtk
from paraview import vtkConstants
from paraview.vtk import *
from paraview.vtk import io

# 
#
#
#
def compute_E_from_L(ncells,grad_v_cell_data):
	#e = numpy.zeros([3,3], float)
	#l = numpy.zeros([3,3], float)
#	ee = [[]*3 for x in xrange(3)]
#	ll = [[]*3 for x in xrange(3)]

	ee = []
	ll = []
	for i in xrange(3):
		ee.append([])
		ll.append([])
		for j in xrange(3):
			ee[i].append(0)
			ll[i].append(0)
	

	#print grad_v_cell_data.GetNumberOfTuples()
	#print grad_v_cell_data.GetNumberOfComponents()


	sr_cell_data = vtk.vtkDoubleArray()
	sr_cell_data.SetName('strain_rate')
	sr_cell_data.SetNumberOfComponents(9)
	sr_cell_data.SetNumberOfTuples(ncells)



	for c in range(0,ncells):
		idx = 9*c

		L = []
		for i in range(0,9):
			L.append(  grad_v_cell_data.GetValue(idx+i) )

		ll[0][0] = L[0]
		ll[0][1] = L[1]
		ll[0][2] = L[2]

		ll[1][0] = L[3]
		ll[1][1] = L[4]
		ll[1][2] = L[5]

		ll[2][0] = L[6]
		ll[2][1] = L[7]
		ll[2][2] = L[8]

		kk = 0
		for i in range(0,3):
			for j in range(0,3):
				ee[i][j] = 0.5 * ( ll[i][j] + ll[j][i] )
				sr_cell_data.SetValue(idx+kk,float(ee[i][j]))
				kk = kk + 1


	return sr_cell_data


def generate_strain_rate(infilename,outfilename):
	# Structured Grid (mesh)
	reader = io.vtkXMLStructuredGridReader()
	reader.SetFileName(infilename)
	reader.Update()
	mesh_in = reader.GetOutput()

	print 'nvertices = ', mesh_in.GetNumberOfPoints()
	print 'ncells    = ', mesh_in.GetNumberOfCells()

	# fetch velocity field
	vel_nodal_data = mesh_in.GetPointData().GetArray('velocity')
	if vel_nodal_data == None:
		print 'Could not locate PointData field named \"velocity\"'

	# fetch gradient field
	grad_v_cell_data = mesh_in.GetCellData().GetArray('VectorGradient')
	if grad_v_cell_data == None:
		print 'Could not locate CellData field named \"VectorGradient\"'


	sr_cell_data = compute_E_from_L(  mesh_in.GetNumberOfCells(), grad_v_cell_data)


	# output mesh
	mesh_in.GetCellData().AddArray( sr_cell_data )

	writer = io.vtkXMLStructuredGridWriter()
	compressor = io.vtkZLibDataCompressor()
	writer.SetCompressor(compressor)

	writer.SetFileName( outfilename )
	writer.SetInput(mesh_in)
	writer.Write()

#	points = mesh_in.GetPoints()
#
#	
#	mesh_out = vtk.vtkStructuredGrid()
#
#	mesh_out.SetNumberOfPoints(  mesh_in.GetNumberOfPoints() )
#	mesh_out.SetNumberOfCells(   mesh_in.GetNumberOfCells() )
#
#	mesh_out.SetPoints( points )
#
#	mesh_out.GetPointData().AddArray( vel_nodal_data )
#	mesh_out.GetCellData().AddArray( sr_cell_data )
#	
#	writer = io.vtkXMLStructuredGridWriter()
#	compressor = io.vtkZLibDataCompressor()
#	writer.SetCompressor(compressor)
#
#	writer.SetFileName( outfilename )
#	writer.SetInput(mesh_out)
#	writer.Write()


def main():
	Connect()

	optparser=OptionParser(usage='usage: %prog -i <filename1>',
												 add_help_option=True,
												 description="""Read vtu file with velocity and grad(V)""" + 
												 """and generate strain rates.""")

	optparser.add_option( "-i", "--input", dest="opt_inputfile",
										help="Input file name", metavar="FILE")

	(options, argv) = optparser.parse_args()

	if options.opt_inputfile == None:
			optparser.print_help()
			sys.exit(1)


	infilename = options.opt_inputfile
	if os.path.splitext(infilename)[1]=='.vts':
		vts_name = os.path.splitext(infilename)[0] + "_sr.vts"
	
		print 'Reading: ' + infilename
		print 'Writing: ' + vts_name
		generate_strain_rate(infilename,vts_name)

if __name__ == '__main__':
	main()
