
import json
import ccodegen2 as C


line = C.include("stdio.h",True)
print(str(line))

func = C.function("myModelCreate",['double a[]','int b'],'int')
lineF = func.open()
lineF += func.statement(['int a = 5;','int b = 44;'])

lineF += func.forBegin('ii',[0,10])
lineF += func.statement(['double valeA = 0.0;'])
lineF += func.forBegin('jj',[0,4])
lineF += func.statement(['double vale = 0.0;'])
lineF += func.statement(['vale = (double)(ii+jj*4);'])
lineF += func.forEnd()
lineF += func.statement(['valeA += (double)(ii*4);'])
lineF += func.forEnd()


lineF += func.ifBegin('i == 6')
lineF += func.statement(['printf(\"i is six\\n\");'])
lineF += func.elseIfBegin('i == 10')
lineF += func.statement(['printf(\"i is ten\\n\");'])
lineF += func.statement(['double vale5[40] = 0.0;'])
lineF += func.forBegin('kk',[0,40])
lineF += func.statement(['vale5[kk] = (double)(ii+jj*4);'])
lineF += func.forEnd()
lineF += func.elseBegin()
lineF += func.statement(['printf(\"i is unknown\\n\");'])
lineF += func.ifEnd()


lineF += func.close()

print(lineF)
print('[Buffered]\n'+func.buffer)

print(func.prototype())




structName = "MPntData"
fname = "cstruct.json"
file = open(fname,'r')
jdata = json.load(file)
file.close()

print(jdata)

for p in jdata['params']:
  print(p)
  print('found: ' + p['name'])
  print('found: ' + p['type'])
  if 'len' not in p:
    print('failed found to find len')

setterList = []
for p in jdata['params']:
  funcName = structName + 'SetField_' + p['name']
  func = C.function(funcName,[structName + ' *p',p['type'] + ' value'],'void')
  func.open()
  func.statement(['p->'+p['name']+' = value;'])
  func.close()
  setterList.append(func)


print(setterList[0].prototype())
print(setterList[1].prototype())
print(setterList[0].write())
print(setterList[1].write())
#print(setter((0)))

