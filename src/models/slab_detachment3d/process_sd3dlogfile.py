

import sys
import re


if len(sys.argv) == 1:
    exit('Error: Requires the full path to the sd3d.logfile to process')

print '# === slab_detachment3d post processing ', sys.argv[1], '=== #'


filename = sys.argv[1]
file = open( filename, 'r' )

outfile = open('sd3d_bmark.dat','w')

#        0       1         2       3     4          5          6          7          8           9           10      11      12      13     14
key = { 'step', 'timeNd', 'time', 'Ys', 'Xn_edge', 'Yn_edge', 'Xn_back', 'Yn_back', 'Zn_front', 'Yn_front', 'Ymin', 'Ymax', 'Vrms', 'Phi', 'vol' }


W0x = 40.0 * 1.0e3 # half slab width
W0z = 250.0 * 1.0e3 # slab length (in the third dimension(
tc = 22.56 * 1.0e6 * ( 3600.0 * 24.0 * 365.0 ); # 22.56 Ma in sec
Vc = W0x / tc
Ly = 660.0 * 1.0e3

# initial volume
V0 = 660.0e3 * 500.0e3 * 500.0e3

scan_remainder = False

for line in file:
    x = line.find('|             step |')
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

        # edge values
        idx = 4 # W - I measured max x coordinate of the slab, so this needs to be rescaled
        W_star = (500.0e3 - float(entries[idx]))/W0x

        idx = 5 # Dn
        DN_star = (Ly - float(entries[idx]))/Ly

        # back face values (y=0 plane in benchmark, z=0 plane in ptatin)
        idx = 6 # W - I measured max x coordinate of the slab, so this needs to be rescaled
        Wn_back_star = (500.0e3 - float(entries[idx]))/W0x

        idx = 7
        DN_back_star = (Ly - float(entries[idx]))/Ly

        # front face values (x=Lx plane in benchmark, z=500 km plane in ptatin)
        idx = 8
        Wn_front_star = (float(entries[idx]))/W0z
        idx = 9
        DN_front_star = (Ly - float(entries[idx]))/Ly

        idx = 10 # Ymin/Ymax
        Zmin_star = (Ly - float(entries[idx]))/Ly
        idx = 11
        Zmax_star = (Ly - float(entries[idx]))/Ly

        idx = 12 # vrms
        Vrms_star = float(entries[idx])/Vc

        idx = 13 # dissipation - SI units
        phi_star = float(entries[idx])

        idx = 14 # vol
        vol_star = (float(entries[idx]) - V0)/V0

#        print '%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e' % (time_star,Ds_star, W_star,DN_star, Wn_back_star,DN_back_star, Wn_front_star,DN_front_star, Zmin_star,Zmax_star, Vrms_star,phi_star,vol_star)

        outfile.write('%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n' % (time_star,Ds_star, W_star,DN_star, Wn_back_star,DN_back_star, Wn_front_star,DN_front_star, Zmin_star,Zmax_star, Vrms_star,phi_star,vol_star))


file.close()
outfile.close()
