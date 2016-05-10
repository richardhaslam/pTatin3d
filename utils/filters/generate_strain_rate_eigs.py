
import math
import sys
from optparse import OptionParser

from paraview.servermanager import *

from paraview import vtk
from paraview import vtkConstants
from paraview.vtk import *
from paraview.vtk import io

#import numpy as NP
#from scipy import *

def CubicSolver_OnlyRealRoots(a2,a1,a0):
	rootList = [0]*3

	Q = (3.0 * a1 - a2 * a2)/9.0
	R = (9.0*a2*a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2)/54.0

	D = Q*Q*Q + R*R
		
	if math.fabs(D) < 1.0e-8:

		if R < 0.0 :
			halfB = - pow( -R, 1.0/3.0 )
		else :
			halfB = pow( R, 1.0/3.0 )

		rootList[0] = -a2/3.0 + 2.0 * halfB
		rootList[1] = -a2/3.0 - halfB
		rootList[2] = -a2/3.0 - halfB
		return rootList
	
	if D > 0.0:
			print "Polynomial discrimanent is positive which means there are complex solutions.\nCannot solve equation"


	theta = math.acos( R/math.sqrt( -Q*Q*Q ) )/3.0
	factor = 2.0 * math.sqrt( -Q )

	rootList[0] = factor * math.cos( theta ) - a2/3.0
	rootList[1] = factor * math.cos( theta + 2.0*math.pi/3.0 ) - a2/3.0
	rootList[2] = factor * math.cos( theta + 4.0*math.pi/3.0 ) - a2/3.0

	return rootList


def SymmetricTensor_CalcAllEigenvalues3D(tensor):

	a2 = - tensor[0][0] - tensor[1][1] - tensor[2][2]
	a1 = (tensor[0][0] * tensor[1][1])  +  \
		 (tensor[0][0] * tensor[2][2])  +  \
		 (tensor[1][1] * tensor[2][2])  -  \
		 (tensor[0][1] * tensor[0][1])  -  \
		 (tensor[0][2] * tensor[0][2])  -  \
		 (tensor[1][2] * tensor[1][2])
	a0 = - (tensor[0][0] * tensor[1][1] * tensor[2][2])  -  \
				2.0 * (tensor[0][1] * tensor[0][2] * tensor[1][2]) +  \
			  (tensor[0][0] * tensor[1][2] * tensor[1][2])  +  \
			  (tensor[1][1] * tensor[0][2] * tensor[0][2])  +  \
			  (tensor[2][2] * tensor[0][1] * tensor[0][1])
	
	rootList = CubicSolver_OnlyRealRoots( a2, a1, a0 )

	return rootList

def EQL(a,b):
	if math.fabs(a-b)<1.0e-8:
		return True
	else:
		return False
		

def SymmetricTensor_CalcEigenvector3D( tensor, eigenvalue ):
	vector = [0]*3

#	print 'ei',eigenvalue

	A = tensor[0][0] - eigenvalue
	B = tensor[1][1] - eigenvalue
	C = tensor[2][2] - eigenvalue
	d = tensor[0][1]
	e = tensor[0][2]
	f = tensor[1][2]

	print A,d,e
	print d,B,f
	print e,f,C,'\n'

	print B*e - f*d, e
	print f*A-d*e, e
	print d*d-B*A,d
	
	
	
	if ((EQL(B*e, f*d)==False) & (EQL( e, 0.0 )==False)):
		print '1'
		vector[0] = 1.0
		vector[1] = (f*A - d*e)/(B*e - f*d)
		vector[2] = (-A - d * vector[1] ) / e 
	elif ((EQL(f*A, d*e)==False) & (EQL( e, 0.0 )==False)) :
		print '2'
		vector[0] = (B*e - f*d)/(f*A - d*e)
		vector[1] = 1.0
		vector[2] = (-d - A*vector[0])/e
	elif ((EQL(d*d, B*A)==False) & (EQL( d, 0.0 )==False)) :
		print '3'
		vector[0] = (B*e - f*d)/(d*d - B*A)
		vector[1] = (-e - A*vector[0])/d
		vector[2] = 1.0
#	else :
#		return False;
	print vector

	#StGermain_VectorNormalise( eigenvector->vector, 3 );
	mag = math.sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] )
	vector[0] = vector[0] / mag
	vector[1] = vector[1] / mag
	vector[2] = vector[2] / mag

	return vector



def SymmetricTensor_CalcEigenvector3D_DAM( tensor, index, eigenvalue ):
	vector = [0]*3
	rhs = [0]*3
	J = []
	iJ = []
	for i in xrange(3):
		J.append([])
		iJ.append([])
		for j in xrange(3):
			J[i].append(0)
			iJ[i].append(0)


	J[0][0] = tensor[0][0] - eigenvalue
	J[1][1] = tensor[1][1] - eigenvalue
	J[2][2] = tensor[2][2] - eigenvalue
	J[0][1] = tensor[0][1]
	J[0][2] = tensor[0][2]
	J[1][2] = tensor[1][2]

	J[2][1] = J[1][2]
	J[2][0] = J[0][2]
	J[1][0] = J[0][1]

	rhs[index] = 1.0
	for kk in range(0,3):
		J[index][kk] = 0.0
	J[index][index] = 1.0


	t4  = J[2][0] * J[0][1]
	t6  = J[2][0] * J[0][2]
	t8  = J[1][0] * J[0][1]
	t10 = J[1][0] * J[0][2]
	t12 = J[0][0] * J[1][1]
	t14 = J[0][0] * J[1][2]
	t17 = 0.1e1 / (t4 * J[1][2] - t6 * J[1][1] - t8 * J[2][2] + t10 * J[2][1] + t12 * J[2][2] - t14 * J[2][1])

	iJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * t17
	iJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * t17
	iJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * t17
	iJ[1][0] = -(-J[2][0] * J[1][2] + J[1][0] * J[2][2]) * t17
	iJ[1][1] = (-t6 + J[0][0] * J[2][2]) * t17
	iJ[1][2] = -(-t10 + t14) * t17
	iJ[2][0] = (-J[2][0] * J[1][1] + J[1][0] * J[2][1]) * t17
	iJ[2][1] = -(-t4 + J[0][0] * J[2][1]) * t17
	iJ[2][2] = (-t8 + t12) * t17

	for i in xrange(3):
		sum = 0.0
		for j in xrange(3):
			sum = sum + iJ[i][j] * rhs[j]
		vector[i] = sum


	mag = math.sqrt( vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2] )
#	if (eigenvalue < 0.0):
#		mag = -mag
		
	vector[0] = vector[0] / mag
	vector[1] = vector[1] / mag
	vector[2] = vector[2] / mag

	return vector


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


	cdata_eigs = vtk.vtkDoubleArray()
	cdata_eigs.SetName('eig_val')
	cdata_eigs.SetNumberOfComponents(3)
	cdata_eigs.SetNumberOfTuples(ncells)


	cdata_v1 = vtk.vtkDoubleArray()
	cdata_v1.SetName('eig_vec1')
	cdata_v1.SetNumberOfComponents(3)
	cdata_v1.SetNumberOfTuples(ncells)

	cdata_v2 = vtk.vtkDoubleArray()
	cdata_v2.SetName('eig_vec2')
	cdata_v2.SetNumberOfComponents(3)
	cdata_v2.SetNumberOfTuples(ncells)

	cdata_v3 = vtk.vtkDoubleArray()
	cdata_v3.SetName('eig_vec3')
	cdata_v3.SetNumberOfComponents(3)
	cdata_v3.SetNumberOfTuples(ncells)


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

		eigs = SymmetricTensor_CalcAllEigenvalues3D(ee)
#		print 'E=',ee
#		print 'eigs = ',eigs

#		ind = argsort(eigs)
#		print ind

#		teigs = [0]*3
#		teigs[0] = eigs[0]
#		teigs[1] = eigs[1]
#		teigs[2] = eigs[2]
#
#		min_ev = 1.0e32
#		min_kk = -1
#		for kk in range(0,3):
#			if teigs[kk] < min_ev:
#				min_ev = teigs[kk]
#				min_kk = kk
#		eigs[0] = teigs[min_kk]
#		teigs[min_kk] = 1.0e32
#
#		min_ev = 1.0e32
#		for kk in range(0,3):
#			if teigs[kk] < min_ev:
#				min_ev = teigs[kk]
#				min_kk = kk
#		eigs[1] = teigs[min_kk]
#		teigs[min_kk] = 1.0e32
#
#		min_ev = 1.0e32
#		for kk in range(0,3):
#			if teigs[kk] < min_ev:
#				min_ev = teigs[kk]
#				min_kk = kk
#		eigs[2] = teigs[min_kk]
#		teigs[min_kk] = 1.0e32
#		
#		print 'sorted',eigs

		v1 = SymmetricTensor_CalcEigenvector3D_DAM(ee,0,eigs[0])
		v2 = SymmetricTensor_CalcEigenvector3D_DAM(ee,1,eigs[1])
		v3 = SymmetricTensor_CalcEigenvector3D_DAM(ee,2,eigs[2])
		
		for kk in range(0,3):
			v1[kk] = v1[kk] * eigs[0]
			v2[kk] = v2[kk] * eigs[1]
			v3[kk] = v3[kk] * eigs[2]



#		print 'v1=',v1
#		print 'v2=',v2
#		print 'v3=',v3

		for kk in range(0,3):
			cdata_eigs.SetValue(3*c+kk,float(eigs[kk]))
			
			cdata_v1.SetValue(3*c+kk,float(v1[kk]))
			cdata_v2.SetValue(3*c+kk,float(v2[kk]))
			cdata_v3.SetValue(3*c+kk,float(v3[kk]))


	return sr_cell_data, cdata_eigs,cdata_v1,cdata_v2,cdata_v3




def compute_Eeigs_from_L(ncells,grad_v_cell_data):
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
	

	cdata_eigs = vtk.vtkDoubleArray()
	cdata_eigs.SetName('eig_val')
	cdata_eigs.SetNumberOfComponents(3)
	cdata_eigs.SetNumberOfTuples(ncells)


	cdata_v1 = vtk.vtkDoubleArray()
	cdata_v1.SetName('eig_vec1')
	cdata_v1.SetNumberOfComponents(3)
	cdata_v1.SetNumberOfTuples(ncells)

	cdata_v2 = vtk.vtkDoubleArray()
	cdata_v2.SetName('eig_vec2')
	cdata_v2.SetNumberOfComponents(3)
	cdata_v2.SetNumberOfTuples(ncells)

	cdata_v3 = vtk.vtkDoubleArray()
	cdata_v3.SetName('eig_vec3')
	cdata_v3.SetNumberOfComponents(3)
	cdata_v3.SetNumberOfTuples(ncells)


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

		for i in range(0,3):
			for j in range(0,3):
				ee[i][j] = 0.5 * ( ll[i][j] + ll[j][i] )

		eigs = SymmetricTensor_CalcAllEigenvalues3D(ee)

		v1 = SymmetricTensor_CalcEigenvector3D_DAM(ee,0,eigs[0])
		v2 = SymmetricTensor_CalcEigenvector3D_DAM(ee,1,eigs[1])
		v3 = SymmetricTensor_CalcEigenvector3D_DAM(ee,2,eigs[2])

		for kk in range(0,3):
			v1[kk] = v1[kk] * eigs[0]
			v2[kk] = v2[kk] * eigs[1]
			v3[kk] = v3[kk] * eigs[2]

#		print 'eigs = ',eigs
		teigs = [0]*3
		teigs[0] = eigs[0]
		teigs[1] = eigs[1]
		teigs[2] = eigs[2]

		min_ev = 1.0e32
		min_kk = -1
		for kk in range(0,3):
			if teigs[kk] < min_ev:
				min_ev = teigs[kk]
				min_kk = kk
		eigs[0] = teigs[min_kk]
		teigs[min_kk] = 1.0e32
		indx0 = min_kk

		min_ev = 1.0e32
		for kk in range(0,3):
			if teigs[kk] < min_ev:
				min_ev = teigs[kk]
				min_kk = kk
		eigs[1] = teigs[min_kk]
		teigs[min_kk] = 1.0e32
		indx1 = min_kk

		min_ev = 1.0e32
		for kk in range(0,3):
			if teigs[kk] < min_ev:
				min_ev = teigs[kk]
				min_kk = kk
		eigs[2] = teigs[min_kk]
		teigs[min_kk] = 1.0e32
		indx2 = min_kk

#		
#		print 'sorted',eigs
#		print indx0,indx1,indx2

		

		for kk in range(0,3):
			cdata_eigs.SetValue(3*c+kk,float(eigs[kk]))
			
			if indx0 == 0:
				cdata_v1.SetValue(3*c+kk,float(v1[kk]))
			elif indx0 == 1:
				cdata_v1.SetValue(3*c+kk,float(v2[kk]))
			else:
				cdata_v1.SetValue(3*c+kk,float(v3[kk]))

			if indx1 == 0:
				cdata_v2.SetValue(3*c+kk,float(v1[kk]))
			elif indx1 == 1:
				cdata_v2.SetValue(3*c+kk,float(v2[kk]))
			else:
				cdata_v2.SetValue(3*c+kk,float(v3[kk]))

			if indx2 == 0:
				cdata_v3.SetValue(3*c+kk,float(v1[kk]))
			elif indx2 == 1:
				cdata_v3.SetValue(3*c+kk,float(v2[kk]))
			else:
				cdata_v3.SetValue(3*c+kk,float(v3[kk]))
			
		if c%500000 == 0:
			print 'Done ' + str(c) + ' cells of ' + str(ncells)


	return cdata_eigs,cdata_v1,cdata_v2,cdata_v3


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


	#sr_cell_data, cdata_eigs,cdata_v1,cdata_v2,cdata_v3 = compute_E_from_L(  mesh_in.GetNumberOfCells(), grad_v_cell_data)
	cdata_eigs,cdata_v1,cdata_v2,cdata_v3 = compute_Eeigs_from_L( mesh_in.GetNumberOfCells(), grad_v_cell_data)

	# output mesh
	#mesh_in.GetCellData().AddArray( sr_cell_data )

	#mesh_in.GetCellData().AddArray( cdata_eigs )

	mesh_in.GetCellData().AddArray( cdata_v1 )
	#mesh_in.GetCellData().AddArray( cdata_v2 )
	#mesh_in.GetCellData().AddArray( cdata_v3 )

#	reader.UpdatePipelineInformation()
#	print reader.PointResultArrayInfo

	writer = io.vtkXMLStructuredGridWriter()
	compressor = io.vtkZLibDataCompressor()
	writer.SetCompressor(compressor)

	writer.SetFileName( outfilename )
	writer.SetInput(mesh_in)
	writer.Write()



#	coords_data = mesh_in.GetPointData().GetArray('coords')
#	vel_data    = mesh_in.GetPointData().GetArray('velocity')
#	points = mesh_in.GetPoints()
#
#	mesh_out = vtk.vtkStructuredGrid()
#
#	mesh_out.SetPoints(points)
#
#	mesh_out.GetPointData().AddArray( coords_data )
#	mesh_out.GetPointData().AddArray( vel_data )
#	mesh_out.GetCellData().AddArray( cdata_v1 )
#	
#	mesh_out.Update()


#	writer = io.vtkXMLStructuredGridWriter()
#	compressor = io.vtkZLibDataCompressor()
#	writer.SetCompressor(compressor)
##
#	writer.SetFileName( 'o.vts' )
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
		vts_name = os.path.splitext(infilename)[0] + "_sr_eigs.vts"
	
		print 'Reading: ' + infilename
		print 'Writing: ' + vts_name
		generate_strain_rate(infilename,vts_name)

if __name__ == '__main__':
	main()
