

import sys
import re


if len(sys.argv) == 1:
    exit('Error: Requires the full path to the sd3d.logfile to process')

#print(sys.argv[1])
print('# ===================== Slab detachment (3D) post processing ===================== #')
print('#  Diagnostics report required by ')
print('#    "A two- and three-dimensional numerical comparison study of slab detachment"')
print('#    C. Thieulot et al., GCubed, 2014')
print('#\n\n')
print('Processing file: ', sys.argv[1])


filename = sys.argv[1]
file = open( filename, 'r' )

outfilename = 'sd3d_diagnostics.dat'
outfile = open(outfilename,'w')

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
        Zmin_star = (float(entries[idx]) - Ly)/Ly
        idx = 11
        Zmax_star = (float(entries[idx]) - Ly)/Ly

        idx = 12 # vrms
        Vrms_star = float(entries[idx])/Vc

        idx = 13 # dissipation - SI units
        phi_star = float(entries[idx])

        idx = 14 # vol
        vol_star = (float(entries[idx]) - V0)/V0

#        print '%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e' % (time_star,Ds_star, W_star,DN_star, Wn_back_star,DN_back_star, Wn_front_star,DN_front_star, Zmin_star,Zmax_star, Vrms_star,phi_star,vol_star)


        # time Ds* W* Dn*   Wn*_back Dn*_back Wn*_front Dn*_front    Zmin,Zmax   Vrms phi vol*
        # 0    1   2  3     4        5        6         7            8    9      10   11  12
        outfile.write('%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n' % (time_star,Ds_star, W_star,DN_star, Wn_back_star,DN_back_star, Wn_front_star,DN_front_star, Zmin_star,Zmax_star, Vrms_star,phi_star,vol_star))


file.close()
outfile.close()

# Open sd3d_bench.dat and read and write specific columns into single files
print('------ Field: Ds ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_Ds.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[1])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[1])) )
outfile.close()
sdfile.close()

print('------ Field: Wxn,Dn-y0 ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_Wxn_y0.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[4])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[4])) )
outfile.close()
sdfile.close()

outfile = open(outfilename,'r')
sdfile = open('ptat3d_Dn_y0.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[5])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[5])) )
outfile.close()
sdfile.close()

print('------ Field: Wzn,Dn-xLx0 ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_Wzn_xLx.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[6])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[6])) )
outfile.close()
sdfile.close()

outfile = open(outfilename,'r')
sdfile = open('ptat3d_Dn_xLx.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[7])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[7])) )
outfile.close()
sdfile.close()

print('------ Field: vrms ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_vrms.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[10])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[10])) )
outfile.close()
sdfile.close()

print('------ Field: phi ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_phi.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[11])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[11])) )
outfile.close()
sdfile.close()

# free surface only
print('------ Field: Z{min,max} ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_zmin.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[8])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[8])) )
outfile.close()
sdfile.close()

outfile = open(outfilename,'r')
sdfile = open('ptat3d_zmax.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[9])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[9])) )
outfile.close()
sdfile.close()

print('------ Field: vol* ------')
outfile = open(outfilename,'r')
sdfile = open('ptat3d_vol.dat','w')
for line in outfile:
    entries = line.split()
    print(entries[12])
    sdfile.write('%.12e %.12e\n' % (float(entries[0]),float(entries[12])) )
outfile.close()
sdfile.close()


