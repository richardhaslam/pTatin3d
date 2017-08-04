
import sys
import json
import ccodegen2 as C


# <depreciated> pre- ccodegen2
def GenEmptyCStructHeader(gname,jdata):
  ccode = []
  ccode.append('')
  ccode.append('')
  ccode.append('#ifndef __' + gname + '_context_h__')
  ccode.append('#define __' + gname + '_context_h__')
  ccode.append('')
  
  ccode.append('typedef struct {')
  ccode.append('  char codegen_signature[256];')
  ccode.append('} ' + gname + 'Context;')
  ccode.append('')
  ccode.append('#endif')
  return(ccode)

# <depreciated> pre- ccodegen2
def GenCStructHeader(gname,jdata):
  ccode = []
  ccode.append('')
  ccode.append('')
  ccode.append('#ifndef __' + gname + '_context_h__')
  ccode.append('#define __' + gname + '_context_h__')
  ccode.append('')

  ccode.append('typedef struct {')
  ccode.append('  char codegen_signature[256];')
  if jdata != None:
    for p in jdata['params']:
      t = p['type']
      n = p['name']
      if 'len' not in p:
        ccode.append('  ' + t + ' ' + n + ';')
      else:
        l = p['len']
        if l == 1:
          print("<< Error >> It is invalid to specify len = 1. len is reserved only for arrays of length > 1");
          sys.exit(1)
        
        ccode.append('  ' + t + ' ' + n + '[' + str(l) + '];')

  ccode.append('} ' + gname + 'Context;')
  ccode.append('')
  ccode.append('#endif')
  return(ccode)


# <depreciated> pre- ccodegen2
def PetscOpts(gname,jdata):
  ccode = []
  ccode.append('PetscErrorCode ierr;');
  ccode.append('PetscBool found;');
  ccode.append('')
  for p in jdata['params']:
    t = p['type']
    n = p['name']
    opt = '-' + gname + '.' + n
    if 'len' not in p:
      if t == 'double' or t == 'PetscReal':
        ccode.append('found = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,\"'+opt+'\",&c->' + n + ',&found);CHKERRQ(ierr);')
      if t == 'int' or t == 'PetscInt':
        ccode.append('found = PETSC_FALSE; ierr = PetscOptionsGetInt(NULL,NULL,\"'+opt+'\",&c->' + n + ',&found);CHKERRQ(ierr);')
    else:
      if t == 'double' or t == 'PetscReal':
        ccode.append('found = PETSC_FALSE; ierr = PetscOptionsGetRealArray(NULL,NULL,\"'+opt+'\",c->' + n + ',&found);CHKERRQ(ierr);')
      if t == 'int' or t == 'PetscInt':
        ccode.append('found = PETSC_FALSE; ierr = PetscOptionsGetIntArray(NULL,NULL,\"'+opt+'\",c->' + n + ',&found);CHKERRQ(ierr);')
  ccode.append('')

  s = ''
  for c in ccode:
    s = s + c + '\n'
  print(s)

  #
  print('[one line joinery]')
  print('\n'.join(ccode))


def PetscOptsByRegion(jdata,classname,prefix,regionId):
  name = classname + "SetFromOptionsByRegion"
  func = C.function(name,[classname + ' *c,PetscInt regionId'],'PetscErrorCode')
  func.petscopen()
  func.statement(['PetscBool found;'])
  func.statement(['char optname[PETSC_MAX_PATH_LEN];'])
  
  for p in jdata['params']:
    t = p['type']
    n = p['name']
    #opt = '-' + prefix + '.' + n + '.r' + str(regionId)
    opt = '-' + prefix + '.' + n + '.r'


    func.statement([''])

    call = 'ierr = PetscSNPrintf(optname,PETSC_MAX_PATH_LEN-1,\"' + opt + '%D\",regionId);CHKERRQ(ierr);'
    func.statement([call])
    
    call = ''
    if 'len' not in p:
      if t == 'double' or t == 'PetscReal':
        call += 'found = PETSC_FALSE; ierr = PetscOptionsGetReal(NULL,NULL,optname,&c->' + n + ',&found);CHKERRQ(ierr);'
      if t == 'int' or t == 'PetscInt':
        call += 'found = PETSC_FALSE; ierr = PetscOptionsGetInt(NULL,NULL,optname,&c->' + n + ',&found);CHKERRQ(ierr);'
    else:
      if t == 'double' or t == 'PetscReal':
        call += 'found = PETSC_FALSE; ierr = PetscOptionsGetRealArray(NULL,NULL,optname,c->' + n + ',&found);CHKERRQ(ierr);'
      if t == 'int' or t == 'PetscInt':
        call += 'found = PETSC_FALSE; ierr = PetscOptionsGetIntArray(NULL,NULL,optname,c->' + n + ',&found);CHKERRQ(ierr);'
    
    func.statement([call])

  func.petscclose()
  return(func.buffer)


def BuildCStruct(jdata,classname):
  struct = C.struct(classname,False)

  ccode = ''
  ccode += struct.open()
  if jdata != None:
    for p in jdata['params']:
      t = p['type']
      n = p['name']
      if 'len' not in p:
        ccode += struct.statement(t + ' ' + n + ';')
      else:
        l = p['len']
        if l == 1:
          print("<< Error >> It is invalid to specify len = 1. len is reserved only for arrays of length > 1");
          sys.exit(1)
        
        ccode += struct.statement(t + ' ' + n + '[' + str(l) + '];')

  ccode += struct.close()
  return(ccode)


fname = "cstruct.json"
file = open(fname,'r')
data = json.load(file)
file.close()

print('[Pickeled JSON data]')
print(data)

if 'className' not in data:
  print('  [Error] Missing essential keyword: \"' + 'className' + '\"')
  sys.exit(1)
if 'classNameShort' not in data:
  print('  [Error] Missing essential keyword: \"' + 'classNameShort' + '\"')
  sys.exit(1)
if 'optionPrefix' not in data:
  print('  [Error] Missing essential keyword: \"' + 'optionPrefix' + '\"')
  sys.exit(1)
if 'params' not in data:
  print('  [Error] Missing essential keyword: \"' + 'params' + '\"')
  sys.exit(1)



print('[Formatted JSON data]')
cname = data['className']
cname_short = data['classNameShort']
optname = data['optionPrefix']

for p in data['params']:
  print(p)
  print('  found: ' + p['name'])
  print('  found: ' + p['type'])
  if 'len' not in p:
    print('  failed found to find len')

#cstructh = GenCStructHeader(cname,data)
#s = ''
#for c in cstructh:
#  s = s + c + '\n'
#print(s)

#PetscOpts(optname,data)


struct = BuildCStruct(data,cname)
func = PetscOptsByRegion(data,cname,optname,2)
print(struct)
print(func)
