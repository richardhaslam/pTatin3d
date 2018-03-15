
import os
import sys
from optparse import OptionParser
import ctypes as ctypes
import struct as struct
import json as json

import databucket_load as db


def main():
  optparser=OptionParser(usage='usage: %prog -i <json_file>',
                         add_help_option=True,
                         description="""Read JSON metadata for """ +
                         """the DataBucket object.""")
    
  optparser.add_option( "-i", "--input", dest="opt_inputfile",
                       help="Input file name", metavar="FILE")

  (options, argv) = optparser.parse_args()
        
  if options.opt_inputfile == None:
    optparser.print_help()
    sys.exit(1)
      

  infilename = options.opt_inputfile
        
  if os.path.splitext(infilename)[1]=='.json':

    print('Reading: ' + infilename)
  
    databucket = db.loadDataBucket(infilename)

    print(databucket)
    print(databucket["MPntStd"].getValues(1))
    print(databucket["MPntPStokes"].getValues(1))
    print(databucket["MPntPStokesPl"].getValues(1))
    try:
      print(databucket["MPntPEnergy"].getValues(1))
    except:
      print("MPntPEnergy was not found")

  else:
    print('Warning: Input file specified is not a JSON file')
    print('Warning: Expected something with the extension .json')
    print('Warning: Found input file: ' + infilename)
    optparser.print_help()
    sys.exit(1)


if __name__ == '__main__':
  main()
