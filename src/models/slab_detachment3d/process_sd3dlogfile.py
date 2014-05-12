

import sys
import re


if len(sys.argv) == 1:
    exit('Error: Requires the full path to the sd3d.logfile to process')

print '# === slab_detachment3d post processing ', sys.argv[1], '=== #'


filename = sys.argv[1]
file = open( filename, 'r' )

#        0       1         2       3     4     5     6       7       8       9      10
key = { 'step', 'timeNd', 'time', 'Ys', 'Yn', 'Xm', 'Ymin', 'Ymax', 'Vrms', 'Phi', 'vol' }


W0 = 40.0 * 1.0e3 # half slab width
tc = 22.56 * 1.0e6 * ( 3600.0 * 24.0 * 365.0 ); # 22.56 Ma in sec
Vc = W0 / tc
Ly = 660.0 * 1.0e3

# initial volume
V0 = 660.0e3 * 500.0e3 * 500.0e3

scan_remainder = False

for line in file:
    x = line.find('# step |')
    if x != -1:
        scan_remainder = True
        continue
    if scan_remainder == True:
        #print line
        entries = line.split()

        idx = 2
        time_star = float(entries[idx])/tc

        idx = 3 # Ds
        Ds_star = (Ly - float(entries[idx]))/Ly

        idx = 4 # Dn
        DN_star = (Ly - float(entries[idx]))/Ly

        idx = 5 # W - I measured max x coordinate of the slab, so this needs to be rescaled
        W_star = (500.0e3 - float(entries[idx]))/W0

        idx = 6 # Ymin/Ymax
        Zmin_star = (Ly - float(entries[idx]))/Ly
        idx = 7
        Z_max = (Ly - float(entries[idx]))/Ly

        idx = 8 # vrms
        Vrms_star = float(entries[idx])/Vc

        idx = 9 # dissipation - scale by what??
        phi_star = float(entries[idx])

        idx = 10 # vol
        vol_star = (float(entries[idx]) - V0)/V0



file.close()