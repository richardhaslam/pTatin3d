
import sys
from math import log as ln

def computeRates(nmeshes,h,E,errName):
  print('-----------------------------------')
  print('         h            ' + errName + '     rate')
  print('%1.4e  %1.10e   ----' % (h[0], E[0]))
  for mesh in range(1,nmeshes):
    Ei   = E[mesh]
    Eim1 = E[mesh-1]
    rate = ln(Ei/Eim1)/ln(h[mesh]/h[mesh-1])
    print('%1.4e  %1.10e  %+.2f' % (h[mesh], Ei, rate))

def computeRatesSummary(nmeshes,h,E,errName):
  rates = []
  rates.append('-')
  for mesh in range(1,nmeshes):
    Ei   = E[mesh]
    Eim1 = E[mesh-1]
    rate = ln(Ei/Eim1)/ln(h[mesh]/h[mesh-1])

    val = '%+.2f' % (rate)

    rates.append(' ' + str(val))
  print('estimated rate ' + errName ,rates)


def parseOutput(filename):
  file = open(filename,'r')
  norms = ''
  for line in file:
    x = line.find('[inorms] ')
    if x != -1:
      norms = line
      break
  file.close()
  return norms


nfiles = len(sys.argv) - 1
if nfiles <= 1:
  print('Error: Require at least 2 input files to estimate a convergence rate')
  exit(1)

raw = []
for f in range(0,nfiles):
  filename = sys.argv[f+1]
  values = parseOutput(filename)
  raw.append(values)

h = []
euL2 = []
euH1s = []
euH1 = []
epL2 = []

for f in range(0,nfiles):
  line = raw[f].lstrip('[inorms] ')
  vals = line.split(',')
  h.append(float(vals[0]));
  euL2.append(float(vals[1]));
  euH1s.append(float(vals[2]));
  euH1.append(float(vals[3]));
  epL2.append(float(vals[4]));

errorNames = [ 'uL2 ' , 'uH1s' , 'uH1 ' , 'pL2 ' ]
errorList = [ euL2 , euH1s , euH1, epL2 ]

for n in range(0,4):
  computeRatesSummary(nfiles,h,errorList[n],errorNames[n])


