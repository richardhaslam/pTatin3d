

import sys
import re


if len(sys.argv) == 1:
	exit('Error: Requires filename to filter')

print '#=== Filtering ', sys.argv[1], '==='

filename = sys.argv[1]
file = open( filename, 'r' )

print '# KSP Component U,V,W,P'
print 'FC_res_a=array(['
fc_cnt = 0
vcycle_cnt = 0
for line in file:

	x = line.find('KSP Component U,V,W,P')
	if x != -1:
		newstr = filter(lambda c: c not in ",", line)
		a = newstr.split()
		#print newstr

		#0 KSP Component U,V,W,P residual norm [ 7.941862009998e+02, 1.234959022576e+01, 3.006820484914e+00, 1.353205878850e-01 ]
		#print a[0], a[7], a[8], a[9], a[10]

		print '[' + a[0] + ',' + str(vcycle_cnt) + ',' + a[7] + ',' + a[8] + ',' + a[9] + ',' + a[10] + '],'

		fc_cnt = fc_cnt + 1


	x = line.find('KSP Residual norm')
	if x != -1:
		newstr = filter(lambda c: c not in ",", line)
		a = newstr.split()
		#print newstr

    #4 KSP Residual norm 2.637068195634e-01
		#print a[0], a[4]

		#print '[' + str(vcycle_cnt) + ',' + a[0] + ',' + a[4] + '],'


		vcycle_cnt = vcycle_cnt + 1


file.close()
print '])'

print '\n'
print 'max_FC_iterations_with_KSP = ' + str(fc_cnt-1)
print 'FCiteration_a = [0]*max_FC_iterations_with_KSP'

print '\n'


####

file = open( filename, 'r' )

fc_max = fc_cnt
fc_cnt = 0
vcycle_cnt = 0
for line in file:

	x = line.find('KSP Component U,V,W,P')
	if x != -1:
		newstr = filter(lambda c: c not in ",", line)
		a = newstr.split()
		#print newstr


		if fc_cnt == (fc_max-1):
			print '])\n'
			break

		if fc_cnt == 0:
			print 'FCiteration_a['+str(fc_cnt)+']=array(['
		else:
			print '])\nFCiteration_a['+str(fc_cnt)+']=array(['


		fc_cnt = fc_cnt + 1


	x = line.find('KSP Residual norm')
	if x != -1:
		newstr = filter(lambda c: c not in ",", line)
		a = newstr.split()
		#print newstr

    #4 KSP Residual norm 2.637068195634e-01
		#print a[0], a[4]

		print '[' + str(vcycle_cnt) + ',' + a[0] + ',' + a[4] + '],'


		vcycle_cnt = vcycle_cnt + 1




file.close()





