#
#  eta = E.eta'
#  x   = L.x'
#  v   = V.v'
#  t   = T.t'
#  p   = P.p'
#  eij = D.eij'
#
#  eta.(d/dx eij ) - dp/dx = rho.g
#  E.(1/L).(V/L) - P/L = rho.g
#              L(u) - grad(p) P/L [ E.(1/L).(V/L) ]^-1 = rho.g [ E.(1/L).(V/L) ]^-1
#  E.(1/L).D  - P/L = rho.g
#
#   
#

import sys
import getopt
#import numpy as np
#import quantities as pq

# Abstract struct class       
class Struct:
	def __init__ (self, *argv, **argd):
		if len(argd):
			# Update by dictionary
			self.__dict__.update (argd)
		else:
			# Update by position
			attrs = filter (lambda x: x[0:2] != "__", dir(self))
			for n in range(len(argv)):
					setattr(self, attrs[n], argv[n])


class PDEStokesScaling (Struct):
	E = 1.0 # viscosity, eta
	L = 1.0 # length
	V = 1.0 # velocity
	T = 1.0 # time
	P = 1.0 # pressure
	D = 1.0 # strain-rate
	rhs_factor = 1.0


class PDEStokesUnits (Struct):
	viscosity = 'Pa.s'
	length = 'm'
	velocity = 'm/s'
	time = 's'
	stress = 'Pa'
	strainrate = '1/s'


class UnitConversions (Struct):
	m2mm = 1000.0
	m2cm = 100.0
	m2km = 1.0/1000.0
	km2m = 1000.0
	s2yr = 1.0/(60.0*60*24.0*365.0)
	yr2s = 60.0*60*24.0*365.0
	m_per_s2cm_per_yr = m2cm / s2yr
	cm_per_yr2m_pers = 1.0 / m_per_s2cm_per_yr


def StokeScalingView(scales,name):
	print '<StokesScales: ' + name + '>'
	print '  eta         ', scales.E
	print '  length      ', scales.L
	print '  velocity    ', scales.V
	print '  time        ', scales.T
	print '  strain-rate ', scales.D
	print '  pressure    ', scales.P
	print '  rhs-scale   ', scales.rhs_factor
	print '</StokesScales>'


def StokeScalingReport_Dimensional2NonDimensional(scales):
	print '<StokesScales: Convert user units into scaled units>'
	print '  eta         -> eta         x ', 1.0/scales.E
	print '  length      -> length      x ', 1.0/scales.L
	print '  velocity    -> velocity    x ', 1.0/scales.V
	print '  time        -> time        x ', 1.0/scales.T
	print '  strain-rate -> strain-rate x ', 1.0/scales.D
	print '  pressure    -> pressure    x ', 1.0/scales.P
	print '</StokesScales>'


# Scaling functions

# compute time, strain-rate, sigma/P, alpha,
def StokesScalingCreate_ELV( _eta, _L, _V ):

	scales  = PDEStokesScaling( E=float(_eta), L=float(_L), V=float(_V) )

	scales.T = scales.L / scales.V

	scales.D = scales.V / scales.L

	fac = scales.E * (1.0/scales.L) * (scales.V/scales.L)
	scales.rhs_factor = 1.0/fac

	scales.P = 1.0 / (fac / scales.L)

	return scales



def ConvertVelocityToSI(value,unit):

	sec_per_year = 60.0 * 60.0 * 24.0 * 365.0
	fac = 1.0

	if   unit == 'm/s':
		fac = 1.0
	elif unit == 'mm/yr':
		fac = (1.0/sec_per_year) * (1.0/1000.0)
	elif unit == 'cm/yr':
		fac = (1.0/sec_per_year) * (1.0/100.0)
	else:
		print 'unknown unit used for velocity'

	sivalue = fac * value
	return sivalue

def ConvertLenthToSI(value,unit):

	fac = 1.0

	if   unit == 'm':
		fac = 1.0
	elif unit == 'km':
		fac = 1000.0
	elif unit == 'cm':
		fac = (1.0/100.0)
	elif unit == 'mm':
		fac = (1.0/1000.0)
	else:
		print 'unknown unit used for length'

	sivalue = fac * value
	return sivalue


def StokesScalingCreate_ConvertUnits(units_from,units_to):

	scales  = PDEStokesScaling()



	return scales


def test_ELV():
	units = UnitConversions()
#	scaling  = PDEStokesScaling(E=2)
#	print scaling.E

#	Vs = ConvertVelocityToSI(3.0,'cm/yr')
#	print '3 cm/yr =>', Vs, 'm/s'
#	print '3 cm/yr =>', 3.0*units.cm_per_yr2m_pers, 'm/s'

	eta_scale      = 400.0
	length_scale   = 100.0
	velocity_scale = 1.0

#	scaling = StokesScalingCreate_ELV( 1.0e24, 1000.0 * units.km2m, 3.0*units.cm_per_yr2m_pers )
	scaling = StokesScalingCreate_ELV( eta_scale, length_scale, velocity_scale )
	StokeScalingView(scaling,'ND')
	StokeScalingReport_Dimensional2NonDimensional(scaling)


def usage_ScalingELV():
	print '[main_ScaleUsingELV]: Provides scaling for Stokes eqation using a viscosity, length and velocity scale as input'
	print 'Usage:'
	print '  -e 4.0 (or --elv_viscosity=4.0) to define viscosity scale'
	print '  -l 4.0 (or --elv_length=4.0) to define length scale'
	print '  -v 4.0 (or --elv_velocity=4.0) to define velocity scale'


#	input_unit_specification = PDEStokesUnits()
#	user_unit_specification  = PDEStokesUnits(length='km')
	# user now modifies all lengths in input to be in km 
	# ...
  #

	# choose the scaling in the units defined by: user_unit_specification 
#	user_output_unit_specification  = PDEStokesUnits(length='km',stress='MPa')

def main_ScaleUsingELV():
	units = UnitConversions()

	# defaults
	eta_scale      = 1.0
	length_scale   = 1.0
	velocity_scale = 1.0

	# parse options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h:e:l:v:", ["help", "elv_viscosity=", "elv_length=", "elv_velocity="])
	except getopt.GetoptError as err:
		print 'Error: Unknown option provided'
		usage_ScalingELV()
		sys.exit(2)

	for o, a in opts:
		if o in   ("-e", "--elv_viscosity"):
			print '  option: --elv_eta = ', a
			eta_scale = a

		elif o in ("-l", "--elv_length"):
			print '  option: --elv_length = ', a
			length_scale = a
			
		elif o in ("-v", "--elv_velocity"):
			print '  option: --elv_velocity = ', a
			velocity_scale = a

		else :
			usage_ScalingELV()
			assert False, "unhandled option"
		
	scaling = StokesScalingCreate_ELV( eta_scale, length_scale, velocity_scale )
	StokeScalingView(scaling,'ND')
	StokeScalingReport_Dimensional2NonDimensional(scaling)


def main1():

	# parse command line options
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
	except getopt.error, msg:
		print msg
		print "for help use --help"
		sys.exit(2)

	# process options
	for o, a in opts:
		if o in ("-h", "--help"):
			print __doc__
			sys.exit(0)
			
	# process arguments
	for arg in args:
		process(arg) # process() is defined elsewhere



def main():

	# test_ELV()

	scale_type = 'ELV'

	if scale_type == 'ELV':
		main_ScaleUsingELV()
	else:
		print 'Unsupported scaling type'
		sys.exit(0)





if __name__ == "__main__":
	main()


