
import sys
import os.path

field_to_fetch = 'Material_Index'

nfiles = 64
for files in range(0,nfiles):
	VTU_FILE_TO_CONVERT = 'picIntegrationPoints.'+str(files)+'.00000.vtu'
	cmd = '/Applications/ParaView3.12.0-RC2.app/Contents/bin/pvpython extract_uw_markers_vtu2binary.py -i ' + str(VTU_FILE_TO_CONVERT) + ' -f ' + str(field_to_fetch)

	os.system( cmd )