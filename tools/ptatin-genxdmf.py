# 
# XDMF PV 4.3.1 Notes
#
# [1] Defining a temporal collection appears uncessary to have a time series displayed in PV 
# [2] Combining many grids into a spatial collection is impractical, users have only access to blocks per time step
# [3] Result of [1] and [2] is that if a temporal collection is desired, different meshes have to be inserted in different root xmf files
# 

import re

def pTatinDMDAStokes_WriteXDMF(mx,my,mz,prefix,time, length_scale, time_scale, velocity_scale, stress_scale):

	nx = 2 * mx + 1
	ny = 2 * my + 1
	nz = 2 * mz + 1

	nodesize = str(nz) + ' ' + str(ny) + ' ' + str(nx)
	cellsize = str(mz) + ' ' + str(my) + ' ' + str(mx)

	fname = prefix+'-StokesGrid.xmf'
	print('  Writing',fname)
	f = open(fname,'w')

	f.write('    <Grid Name="Stokes" GridType="Uniform">\n')

	time_value = float(time) * float(time_scale['Scale'])
	f.write('      <Time Value="' + str(time_value) + '"/>\n')
	f.write('      <Topology TopologyType="3DSMesh" Dimensions="' + nodesize + '"/>\n')
	f.write('      <Geometry GeometryType="XYZ">\n')
	scale = length_scale['Scale']
	f.write('        <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + nodesize + ' 3">\n')
	f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + nodesize + ' 3">\n')
	f.write('            ' + prefix + '.dmda-velocity.coords\n')
	f.write('          </DataItem>\n')
	f.write('        </DataItem>\n')
	f.write('      </Geometry>\n')
	f.write('\n')

	fieldname = 'velocity' + ' [' + velocity_scale['Unit'] + ']'
	scale = velocity_scale['Scale']
	f.write('      <Attribute Name="'+ fieldname  +'" AttributeType="Vector" Center="Node">\n')
	f.write('        <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + nodesize + ' 3">\n')
	f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + nodesize + ' 3">\n')
	f.write('            ' + prefix + '.dmda-Xu\n')
	f.write('          </DataItem>\n')
	f.write('        </DataItem>\n')
	f.write('      </Attribute>\n')
	f.write('    </Grid>\n')
	f.write('\n')

	f.close()

def pTatinDMDAStokesP0_WriteXDMF(mx,my,mz,prefix,time, length_scale, time_scale, stress_scale):

	q2nx = 2*mx + 1
	q2ny = 2*my + 1
	q2nz = 2*mz + 1
	q2nodesize = str(q2nz) + ' ' + str(q2ny) + ' ' + str(q2nx)

	nx = mx + 1
	ny = my + 1
	nz = mz + 1

	nodesize = str(nz) + ' ' + str(ny) + ' ' + str(nx)
	cellsize = str(mz) + ' ' + str(my) + ' ' + str(mx)

	fname = prefix+'-StokesP0Grid.xmf'
	print('  Writing',fname)
	f = open(fname,'w')

	f.write('<Grid Name="StokesP0_' + prefix + '" GridType="Uniform">\n')

	time_value = float(time) * float(time_scale['Scale'])
	f.write('  <Time Value="' + str(time_value) + '"/>\n')
	f.write('\n')

	f.write('  <Topology TopologyType="3DSMesh" Dimensions="' + nodesize + '"/>\n')
	f.write('  <Geometry GeometryType="XYZ">\n')

	scale = length_scale['Scale']
	f.write('    <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + nodesize + ' 3">\n')
	f.write('      <DataItem Reference="XML">\n')
	f.write('      /Xdmf/Domain/Grid[@Name="TimeSeries-StokesP0Grid"]/Grid[@Name="StokesP0_' + prefix + '"]/DataItem[@Name="subset_coor"]\n')
	f.write('      </DataItem>\n')
	f.write('    </DataItem>\n')

	f.write('  </Geometry>\n')
	f.write('\n')


	fieldname = 'pressure' + ' [' + stress_scale['Unit'] + ']'
	scale = stress_scale['Scale']
	f.write('  <Attribute Name="' + fieldname + '" AttributeType="Scalar" Center="Cell">\n')
	f.write('    <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + cellsize + '">\n')
	f.write('      <DataItem Reference="XML">\n')
	f.write('      /Xdmf/Domain/Grid[@Name="TimeSeries-StokesP0Grid"]/Grid[@Name="StokesP0_' + prefix + '"]/DataItem[@Name="subset_pressure"]\n')
	f.write('      </DataItem>\n')
	f.write('    </DataItem>\n')
	f.write('  </Attribute>\n')
	f.write('\n')

	f.write('  <!-- extract coords -->\n')
	f.write('  <DataItem Name="subset_coor" ItemType="HyperSlab" Dimensions="' + nodesize + ' 3" Type="HyperSlab">\n')
	f.write('    <DataItem Dimensions="3 4" Format="XML">\n')
	f.write('      0 0 0 0\n')
	f.write('      2 2 2 1\n')
	f.write('      ' + nodesize + ' 3\n')
	f.write('    </DataItem>\n')
	f.write('    <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + q2nodesize + ' 3">\n')
	f.write('      ' + prefix + '.dmda-velocity.coords\n')
	f.write('    </DataItem>\n')
	f.write('  </DataItem>\n')
	f.write('\n')

	f.write('  <!-- extract pressure -->\n')
	f.write('  <DataItem Name="subset_pressure" ItemType="HyperSlab" Dimensions="' + cellsize + '" Type="HyperSlab">\n')
	f.write('    <DataItem Dimensions="3 4" Format="XML">\n')
	f.write('      0 0 0 0\n')
	f.write('      1 1 1 4\n')
	f.write('      ' + cellsize + ' 1\n')
	f.write('    </DataItem>\n')
	f.write('    <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + cellsize + ' 4">\n')
	f.write('      ' + prefix + '.dmda-Xp\n')
	f.write('    </DataItem>\n')
	f.write('  </DataItem>\n')
	f.write('\n')

	f.write('</Grid>\n')
	f.close()

def pTatinDMDAEnergy_WriteXDMF(mx,my,mz,prefix,time, length_scale, time_scale, temp_scale):

	nx = mx + 1
	ny = my + 1
	nz = mz + 1

	nodesize = str(nz) + ' ' + str(ny) + ' ' + str(nx)

	fname = prefix+'-EnergyGrid.xmf'
	print('  Writing',fname)
	f = open(fname,'w')

	f.write('    <Grid Name="Energy" GridType="Uniform">\n')

	time_value = float(time) * float(time_scale['Scale'])
	f.write('      <Time Value="' + str(time_value) + '"/>\n')
	f.write('      <Topology TopologyType="3DSMesh" Dimensions="' + nodesize + '"/>\n')
	f.write('      <Geometry GeometryType="XYZ">\n')

	scale = length_scale['Scale']
	f.write('        <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + nodesize + ' 3">\n')
	f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + nodesize + ' 3">\n')
	f.write('            ' + prefix + '.dmda-energy.coords.vec\n')
	f.write('          </DataItem>\n')
	f.write('        </DataItem>\n')
	f.write('      </Geometry>\n')
	f.write('\n')

	fieldname = 'temperature' + ' [' + temp_scale['Unit'] + ']'
	f.write('      <Attribute Name="'+ fieldname + '" AttributeType="Scalar" Center="Node">\n')
	f.write('        <DataItem ItemType="Function" Function="1.0 * $0" Dimensions="' + nodesize + ' 1">\n')
	f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + nodesize + ' 1">\n')
	f.write('            ' + prefix + '.dmda-energy.temperature.vec\n')
	f.write('          </DataItem>\n')
	f.write('        </DataItem>\n')
	f.write('      </Attribute>\n')
	f.write('    </Grid>\n')
	f.close()


def pTatinDMDACell_WriteXDMF(mx,my,mz,prefix,time, length_scale, time_scale, active_variables):

	nx = mx + 1
	ny = my + 1
	nz = mz + 1

	nodesize = str(nz) + ' ' + str(ny) + ' ' + str(nx)
	cellsize = str(mz) + ' ' + str(my) + ' ' + str(mx)

	fname = prefix+'-MPCellGrid.xmf'
	print('  Writing',fname)
	f = open(fname,'w')

	f.write('    <Grid Name="Cell" GridType="Uniform">\n')

	time_value = float(time) * float(time_scale['Scale'])
	f.write('      <Time Value="' + str(time_value) + '"/>\n')
	f.write('      <Topology TopologyType="3DSMesh" Dimensions="' + nodesize + '"/>\n')
	f.write('      <Geometry GeometryType="XYZ">\n')
	
	scale = length_scale['Scale']
	f.write('        <DataItem ItemType="Function" Function="' + str(scale) + ' * $0" Dimensions="' + nodesize + ' 3">\n')
	f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + nodesize + ' 3">\n')
	f.write('            ' + prefix + '.dmda-cell.coords.vec\n')
	f.write('          </DataItem>\n')
	f.write('        </DataItem>\n')
	f.write('      </Geometry>\n')

	list, scale = [],[]

	for var in active_variables:
		if var in ['region','viscosity','density','plastic_strain','yield_indicator','diffusivity','heat_source']:
			if active_variables[var]['Active'] is True:
				list.append('%s [%s]'%(var,active_variables[var]['Unit']))
				scale.append(active_variables[var]['Scale'])
		
	index = 0
	for member in list:
		scale_member = scale[index]
		f.write('\n')
		f.write('      <Attribute Name="' + member + '" AttributeType="Scalar" Center="Cell">\n')
		f.write('        <DataItem ItemType="Function" Function="' + str(scale_member) +' * $0" Dimensions="' + cellsize + ' 1">\n')
		f.write('          <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Big" Seek="8" Dimensions="' + cellsize + ' 1">\n')
		f.write('            ' + prefix + '.dmda-cell.' + member.split()[0]+ '.vec\n')
		f.write('          </DataItem>\n')
		f.write('        </DataItem>\n')
		f.write('      </Attribute>\n')
		index = index + 1

	f.write('    </Grid>\n')
	f.close()


#
#<?xml version="1.0"?>
#<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
#<Collection>
#  <DataSet timestep="0.000000e+00" file="./step000000_X.pvts"/>
#  <DataSet timestep="1.000000e-13" file="./step000001_X.pvts"/>
#
def pTatinXDMF_FilterPVDContents(filename):
	print('Filter temporal data ' + filename )
	# for each time/name pair
	#	write three xdmf files
	#	include names in Temporal collection
	
	timeseries = []

	file = open(filename,'r')

	# skip three lines
	file.readline()
	file.readline()
	file.readline()

	for line in file:

		res = line.split('"')[1::2]
		if not (res):
			continue

		time_value = res[0]
		step_value = -1		

		matchObj = re.match( r'(.*)step([0-9]{6})(.*)', res[1] )
		if matchObj:
#			print("m.g()", matchObj.group())
#			print("m.g(2)", matchObj.group(2))
			step_value = matchObj.group(2)
		else:
			print('Failed to match step pattern associated with file name',res[1])
			exit(0)

		print('  time =', time_value, '; step index =', step_value ,'; data file =',res[1])

		timeseries.append([step_value, time_value])
	file.close()

	return(timeseries)


def pTatinXDMF_WritePerStepFiles(timeseries,timeseries_energy,timeseries_mptCell,mx,my,mz,units):
	print('Writing per-step xdmf meta')

	# length_scale, time_scale, velocity_scale, stress_scale)
	ls = units['length']
	ts = units['time']
	vs = units['velocity']
	ss = units['stress']
	temps = units['temperature']

	for time_entry in timeseries:
		pTatinDMDAStokes_WriteXDMF(mx,my,mz,'step'+ time_entry[0], time_entry[1], ls,ts,vs,ss)

	for time_entry in timeseries:
		pTatinDMDAStokesP0_WriteXDMF(mx,my,mz,'step'+ time_entry[0], time_entry[1], ls,ts,ss)

	# length_scale, time_scale
	for time_entry in timeseries_energy:
		pTatinDMDAEnergy_WriteXDMF(mx,my,mz,'step'+ time_entry[0], time_entry[1], ls,ts,temps)

	# length_scale, time_scale, units dict
	for time_entry in timeseries_mptCell:
		pTatinDMDACell_WriteXDMF(mx,my,mz,'step'+ time_entry[0], time_entry[1], ls,ts,units)


def pTatinXDMF_WriteInclude(timeseries, gridname_list):
	print('Writing per-step xdmf')

	for time_entry in timeseries:
		fname = 'ptatin-step' + time_entry[0] + '.xmf'
		print('Writing include file',fname)
		f = open(fname,'w')
		for gridname in gridname_list:
			incfile = 'step'+ time_entry[0] + '-' + gridname + '.xmf'
			print('  Include filename', incfile)
			f.write('      <xi:include href=" + incfile + "/>\n')
		f.close()


def pTatinXDMF_WriteIndividualTemporalCollection(timeseries, gridname_list,suffix):
	print('Write temporal collection')

	for gridname in gridname_list:

		if suffix == None:
			fname = 'ptatin-timeseries-' + gridname + '.xmf'
		else:
			fname = 'ptatin-timeseries-' + gridname + suffix + '.xmf'

		print('  Writing',fname)
		f = open(fname,'w')

		f.write('<?xml version="1.0" ?>\n')
		f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
		f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')
		f.write('  <Domain>\n')

		timeseries_name = 'TimeSeries-' + gridname
		#f.write('<!--\n')
		f.write('    <Grid Name="' + timeseries_name + '" GridType="Collection" CollectionType="Temporal">\n')
		#f.write('-->\n')

		for time_entry in timeseries:
			f.write('      <xi:include href="step'+ time_entry[0] + '-' + gridname + '.xmf"/>\n')

		#f.write('<!--\n')
		f.write('    </Grid>\n')
		#f.write('-->\n')

		f.write('  </Domain>\n')
		f.write('</Xdmf>\n')

		f.close()


def pTatinCreateUnitDictionary():

	# Load quantity data from a file generated by ptatin

	L = 1.0e5
	V = 1.0e-10
	E = 1.0e26

	T = L/V
	STR = V/L
	P = E/T
	D = L*L/T
	RHO = -100000.0 #P * T*T/L
	
	quantities = [
	# ------------------------------------
	            { 'Quantity': 'length',
	              'Unit':     'm',
	              'Scale':    '%e'%L, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'velocity',
	              'Unit':     'm/s',
	              'Scale':    '%e'%V, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'viscosity',
	              'Unit':     'Pa.s',
	              'Scale':    '%e'%E, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'time',
	              'Unit':     's',
	              'Scale':    '%e'%T, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'stress',
	              'Unit':     'Pa',
	              'Scale':    '%e'%P, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'density',
	              'Unit':     'kg/m^3',
	              'Scale':    '%e'%RHO, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'strain-rate',
        	      'Unit':     '1/s',
        	      'Scale':    '%e'%STR, 
	              'Active':	   True },
	# ------------------------------------
        	    { 'Quantity': 'diffusivity',
        	      'Unit':     'm^2/s',
        	      'Scale':    '%e'%D, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'temperature',
	              'Unit':     'degC',
	              'Scale':    1.0, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'heat_source',
	              'Unit':     'K/s',
	              'Scale':    '%e'%(1.0/T), 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'heat-flux',
	              'Unit':     'W/m^2',
	              'Scale':    1.0, 
	              'Active':	   True },
	# ------------------------------------
	            { 'Quantity': 'region',
	              'Unit':     '-',
	              'Scale':    1.0, 
	              'Active':	   True },	
	# ------------------------------------
	            { 'Quantity': 'plastic_strain',
	              'Unit':     '-',
	              'Scale':    1.0, 
	              'Active':	   True },
	# ------------------------------------
				{ 'Quantity': 'yield_indicator',
	              'Unit':     '-',
	              'Scale':    1.0, 
	              'Active':	   True },
	]

	# See URL
	#  http://stackoverflow.com/questions/16761884/python-access-list-element-by-key
	# for explaination

	units = {x['Quantity']: x for x in quantities}
	
	#X = units['velocity']
	#print(units['length'])
	print(units)

	return(units)


# GeoUnits defines the following:
#   "stress"   MPa
#   "time"     kyr
#   "velocity" cm/yr
#   "length"   km       
def pTatinCreateGeodynamicUnitDictionary(unit):

	# Convert units into km,cm/kyr,yr
	units_geo = units

	X = units_geo['length']
	X['Scale'] = float(X['Scale']) / 1.0e3
	X['Unit'] = 'km'

	X = units_geo['time']
	sec_p_yr = 60.0 * 60.0 * 24.0 * 365.0
	sec_p_kyr = sec_p_yr * 1000.0 	
	X['Scale'] = float(X['Scale']) / sec_p_kyr
	X['Unit'] = 'kyr'

	X = units_geo['velocity']
	factor = 1.0e-2 / sec_p_yr
	X['Scale'] = float(X['Scale']) / factor
	X['Unit'] = 'cm/yr'

	X = units_geo['stress']
	factor = 1.0e6
	X['Scale'] = float(X['Scale']) / factor
	X['Unit'] = 'MPa'

	#print( '"' + X['Quantity'] + '"' , 'reference value =', X['Scale'] , '[' + X['Unit'] + ']')

	print(units_geo)

	return(units_geo)


# ___________________________________________________________
# ___________________________________________________________
#

# Specify Q2 mesh sizes
mx = 512
my = 128
mz = 512

# Specify which xdmf output to generate from the dmda
Stokes = True
Energy = True
MptCell = True

# Specify which of the following MptCell quantities are available
region = True
viscosity = True
density = True
plastic_strain = True
yield_indicator  = False
diffusivity = True
heat_source = True



units = pTatinCreateUnitDictionary()
geounits = pTatinCreateGeodynamicUnitDictionary(units)


if Stokes is True:
	print '\n\tProcessing Velocity and Pressure'
	timeseries = pTatinXDMF_FilterPVDContents( 'timeseries_X.pvd' )
	gridname_list = [ 'StokesGrid', 'StokesP0Grid' ]
	pTatinXDMF_WritePerStepFiles( timeseries, [], [], mx,my,mz, geounits )
	pTatinXDMF_WriteIndividualTemporalCollection( timeseries, gridname_list, '_geoSI' )

if Energy is True:
	print '\n\tProcessing Temperature'
	timeseries_energy = pTatinXDMF_FilterPVDContents('timeseries_energy.pvd')
	gridname_list_energy = [ 'EnergyGrid' ]
	pTatinXDMF_WritePerStepFiles( [], timeseries_energy, [], mx,my,mz, geounits )
	pTatinXDMF_WriteIndividualTemporalCollection( timeseries_energy, gridname_list_energy, '_geoSI' )

if MptCell is True:
	print '\n\tProcessing Material Point Cell data'
	timeseries_mptCell = pTatinXDMF_FilterPVDContents('timeseries_mpoints_cell.pvd')
	gridname_list_mptCell = [ 'MPCellGrid' ]
	#update the units dictionary to discard undesired variables
	geounits['region']['Active'] = region
	geounits['viscosity']['Active'] = viscosity
	geounits['density']['Active'] = density
	geounits['plastic_strain']['Active'] = plastic_strain
	geounits['yield_indicator']['Active'] = yield_indicator
	geounits['diffusivity']['Active'] = diffusivity
	geounits['heat_source']['Active'] = heat_source
	pTatinXDMF_WritePerStepFiles( [], [], timeseries_mptCell, mx,my,mz, geounits )
	pTatinXDMF_WriteIndividualTemporalCollection( timeseries_mptCell, gridname_list_mptCell, '_geoSI' )


#pTatinXDMF_WriteIndividualTemporalCollection( timeseries, gridname_list, None )





