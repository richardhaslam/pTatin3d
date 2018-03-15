
import os
import sys
from optparse import OptionParser
import ctypes as ctypes
import struct as struct
import json as json


# MPntStd ===================================
def unpack_MPntStd_allAtOnce(f,L):
  bytes_int   = ctypes.sizeof(ctypes.c_int)
  bytes_long   = ctypes.sizeof(ctypes.c_long)
  bytes_double = ctypes.sizeof(ctypes.c_double)
  
  bytes_set = bytes_long + 3*bytes_double + 3*bytes_double + bytes_int + bytes_int
  
  entry_byte = f.read(bytes_set*L)
  
  points = struct.unpack("lddddddii"*L, entry_byte)
  return points

class MPntStd:
  def __init__(self):
    
    bytes_int    = ctypes.sizeof(ctypes.c_int)
    bytes_long   = ctypes.sizeof(ctypes.c_long)
    bytes_double = ctypes.sizeof(ctypes.c_double)
    
    bytes_set = bytes_long + 3*bytes_double + 3*bytes_double + bytes_int + bytes_int
    self.point_bytes = bytes_set
    self.numberFlatMembers = 1 + 3 + 3 + 1 + 1
  
  def load(self,file,L,atomicSize):
    self.L = L
    
    if self.point_bytes != atomicSize:
      print('  ** atomic size mismatch **')
      sys.exit(1)
    
    self.byte_array = unpack_MPntStd_allAtOnce(file,L)
  
  def getValues(self,index):
    s = index * self.numberFlatMembers
    e = s + self.numberFlatMembers
    point = self.byte_array[ s : e ]
    return point
  
  def getMemberCSizes(self):
    sizes = [[1 , "long int"],
             [3 , "double"],
             [3 , "double"],
             [1 , "int"],
             [1 , "int"]]
    return sizes
  
  def getMemberNames(self):
    names = ["point_index",
             "global_coord",
             "local_coord",
             "phase_index",
             "local_element_index"]
    return names

# MPntPStokes ===================================
def unpack_MPntPStokes_allAtOnce(f,L):
  bytes_double = ctypes.sizeof(ctypes.c_double)
  
  bytes_set = 2*bytes_double
  entry_byte = f.read(bytes_set*L)
  points = struct.unpack("dd"*L, entry_byte)
  return points

class MPntPStokes:
  def __init__(self):
    
    bytes_double = ctypes.sizeof(ctypes.c_double)
    
    bytes_set = 2*bytes_double
    self.point_bytes = bytes_set
    self.numberFlatMembers = 1 + 1
  
  def load(self,file,L,atomicSize):
    self.L = L
    
    if self.point_bytes != atomicSize:
      print('  ** atomic size mismatch **')
      sys.exit(1)
    
    self.byte_array = unpack_MPntPStokes_allAtOnce(file,L)
  
  def getValues(self,index):
    s = index * self.numberFlatMembers
    e = s + self.numberFlatMembers
    point = self.byte_array[ s : e ]
    return point
  
  def getMemberCSizes(self):
    sizes = [[1 , "double"],
             [1 , "double"]]
    return sizes
  
  def getMemberNames(self):
    names = ["eta_effective",
             "density"]
    return names

# MPntPStokesPl ===================================
def unpack_MPntPStokesPl_allAtOnce(f,L):
  bytes_float = ctypes.sizeof(ctypes.c_float)
  bytes_short = ctypes.sizeof(ctypes.c_short)
  
  bytes_set = bytes_float + bytes_short
  entry_byte = f.read(bytes_set*L)
  points = struct.unpack("fh"*L, entry_byte)
  return points

# call this method of the short appears to have been byte aligned to 4 bytes
def unpack_MPntPStokesPl_padded_allAtOnce(f,L):
  bytes_float = ctypes.sizeof(ctypes.c_float)
  bytes_short = ctypes.sizeof(ctypes.c_short)
  
  bytes_set = bytes_float + bytes_float
  entry_byte = f.read(bytes_set*L)
  points = struct.unpack("ff"*L, entry_byte)
  return points

class MPntPStokesPl:
  def __init__(self):
    
    bytes_float = ctypes.sizeof(ctypes.c_float)
    bytes_short = ctypes.sizeof(ctypes.c_short)
    
    bytes_set = bytes_float + bytes_short
    self.point_bytes = bytes_set
    self.numberFlatMembers = 1 + 1
    self.aligned = False
  
  def load(self,file,L,atomicSize):
    self.L = L
    
    if self.point_bytes != atomicSize:
      self.aligned = True
      print('  ** atomic size mismatch **')
      print('  ** Expected ',self.point_bytes,' from python def')
      print('  ** Expected ',atomicSize,' from JSON def')
      print('  ** Likely cause is 4-byte alignment when sizeof(struct MPntPStokesPl) was used to define the datafield size')
      print('  ** Loading assuming float,float')
      self.byte_array = unpack_MPntPStokesPl_padded_allAtOnce(file,L)
    else:
      self.byte_array = unpack_MPntPStokesPl_allAtOnce(file,L)

  def getValues(self,index):
    s = index * self.numberFlatMembers
    e = s + self.numberFlatMembers
    point = self.byte_array[ s : e ]
    return point
  
  def getMemberCSizes(self):
    sizes = [[1 , "float"],
             [1 , "short"]]
    if self.aligned == True:
      sizes = [[1 , "float"],
               [1 , "short-->float"]]
         
    return sizes
  
  def getMemberNames(self):
    names = ["plastic_strain",
             "yield_indicator"]
    return names

# MPntPEnergy ===================================
def unpack_MPntPEnergy_allAtOnce(f,L):
  bytes_double = ctypes.sizeof(ctypes.c_double)
  
  bytes_set = 2*bytes_double
  entry_byte = f.read(bytes_set*L)
  points = struct.unpack("dd"*L, entry_byte)
  return points

class MPntPEnergy:
  def __init__(self):
    
    bytes_double = ctypes.sizeof(ctypes.c_double)
    
    bytes_set = 2*bytes_double
    self.point_bytes = bytes_set
    self.numberFlatMembers = 1 + 1
  
  def load(self,file,L,atomicSize):
    self.L = L
    
    if self.point_bytes != atomicSize:
      print('  ** atomic size mismatch **')
      sys.exit(1)
    
    self.byte_array = unpack_MPntPEnergy_allAtOnce(file,L)
  
  def getValues(self,index):
    s = index * self.numberFlatMembers
    e = s + self.numberFlatMembers
    point = self.byte_array[ s : e ]
    return point
  
  def getMemberCSizes(self):
    sizes = [[1 , "double"],
             [1 , "double"]]
    return sizes
  
  def getMemberNames(self):
    names = ["diffusivity",
             "heat_source"]
    return names

# GenericData ===================================
def read_bytes_allAtOnce(f,total_bytes_to_read):
  bytes = f.read(total_bytes_to_read)
  return bytes

def unpack_double_allAtOnce(entry_byte,L,bs):
  format = "d"*bs*L
  points = struct.unpack(format, entry_byte)
  return points

class GenericData:
  def __init__(self,name,totalBytes,L):
    self.name = name
    self.totalBytes = totalBytes
    self.L = L
  
  def load(self,file):
    self.bytes = read_bytes_allAtOnce(file,self.totalBytes)

  def unpack(self,blockSize,cType):
    self.blockSize = blockSize
    self.cType = cType
    if cType == ctypes.c_double:
      self.byte_array = unpack_double_allAtOnce(self.bytes,self.L,self.blockSize)
      self.bytes = []
    else:
      print('*** only cType = double is supported')
      sys.exit(1)

  def getValues(self,index):
    s = index * self.blockSize
    e = s + self.blockSize
    point = self.byte_array[ s : e ]
    return point
  
  def getMemberCSizes(self):
    if self.cType == ctypes.c_double:
      sizes = [[self.blockSize , "double"]]
    else:
      print('*** only cType = double is supported')
      sys.exit(1)
    return sizes


def loadDataBucket(inputfile):

  jf = open(inputfile, "r")
  json_meta = json.load(jf)
  jf.close()

  #print(json_meta["DataBucket"]["partition"]["length"])
  csize = json_meta["DataBucket"]["partition"]["commSize"]
  lengths = json_meta["DataBucket"]["partition"]["length"]
  tlen = 0
  for r in range(0,csize):
    tlen += lengths[r]


  nfields = json_meta["DataBucket"]["nfields"]
  fields = json_meta["DataBucket"]["fields"]
  metaField = []
  for item in range(0,nfields):
    fieldName = fields[item]["fieldName"]
    atomicSize = fields[item]["atomicSize"]
    fieldFileName = fields[item]["fileName"]
    data = { "fieldName": fieldName, "atomicSize": atomicSize, "fileName": fieldFileName }
    metaField.append(data)


  # check all the file names are actually the same
  for f in range(1,nfields):
    if fields[item]["fileName"] != fields[0]["fileName"]:
      print('Reader is only valid if all DataBucket fields are written into a single file')
      sys.exit(1)

  relative_path_filename = fields[item]["fileName"]
  dirs = relative_path_filename.split('/')
  filename = dirs[-1]
  #print('Loading binary file: ',filename)

  dirs = inputfile.split('/')
  filename_p = ""
  ndirs = len(dirs)
  for d in range(0,ndirs-1):
    filename_p += dirs[d] + '/'
  filename_p += filename
  print('Loading binary file: ',filename_p)

  bytes_long   = ctypes.sizeof(ctypes.c_long)

  databucket = {}

  dbfile = open(filename_p, "rb")

  for f in range(0,nfields):
    
    print('Loading field named: ' + metaField[f]["fieldName"])
    header_byte = dbfile.read(bytes_long)
    total = struct.unpack("l", header_byte)[0]
    print('  bytes total to load ', total)
    
    if metaField[f]["fieldName"] == "MPntStd":
      print('  matched MPntStd')
      point = MPntStd()
      point.load(dbfile,tlen,metaField[f]["atomicSize"])
      
      # Merge dict
      loaded = {"MPntStd":point}
      databucket = dict(databucket, **loaded)

    elif metaField[f]["fieldName"] == "MPntPStokes":
      print('  matched MPntPStokes')
      point = MPntPStokes()
      point.load(dbfile,tlen,metaField[f]["atomicSize"])
      
      loaded = {"MPntPStokes":point}
      databucket = dict(databucket, **loaded)

    elif metaField[f]["fieldName"] == "MPntPStokesPl":
      print('  matched MPntPStokesPl')
      point = MPntPStokesPl()
      point.load(dbfile,tlen,metaField[f]["atomicSize"])
      
      loaded = {"MPntPStokesPl":point}
      databucket = dict(databucket, **loaded)
    
    elif metaField[f]["fieldName"] == "MPntPEnergy":
      print('  matched MPntPEnergy')
      point = MPntPEnergy()
      point.load(dbfile,tlen,metaField[f]["atomicSize"])
      
      loaded = {"MPntPEnergy":point}
      databucket = dict(databucket, **loaded)
    
    else:
      print('An unknown data type was encounted. Loading raw bytes')
      point = GenericData(metaField[f]["fieldName"],total,tlen)
      point.load(dbfile)

      loaded = {metaField[f]["fieldName"]:point}
      databucket = dict(databucket, **loaded)

  dbfile.close()

  return databucket
