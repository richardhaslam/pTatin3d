
import os
import shutil as shell
import re

mname = 'SlabEPSL'

# update registration
def ptatinModelRegUpdate():
  shell.copyfile('ptatin_models_reg.c','ptatin_models_reg_dev.c')

  c_prototype = []
  ptatin_reg_list = []
  for dir in os.listdir("sandbox"):
    print("[Found model] " + os.path.join("./sandbox", dir))

    functionName = 'pTatinModelCreate_' + dir
    c_prototype.append('extern PetscErrorCode ' + functionName + '(pTatinModel);')
    
    stringName = dir
    for f in re.findall("([A-Z]+)", stringName):
      stringName = stringName.replace(f, '_'+f.lower())
    stringName = stringName.lstrip('_')
    print(stringName)
  
    
    ptatin_reg_list.append('ierr = pTatinModelDynamicRegister("'+stringName+'",'+functionName+');CHKERRQ(ierr);')

  print(c_prototype)
  print(ptatin_reg_list)

  file = open('ptatin_models_reg_dev.c','r')
  fileN = open('ptatin_models_reg_devM.c','w')
  for line in file:
    flat = line.lstrip()
    flat = flat.rstrip()

    fileN.write(line)
    if flat == '/* sandbox-model-register-prototypes */':
      for l in c_prototype:
        fileN.write(l + '\n')

    if flat == '/* sandbox-model-register */':
      for l in ptatin_reg_list:
        fileN.write('  ' + l + '\n')

  file.close()
  fileN.close()

  shell.move('ptatin_models_reg_devM.c','ptatin_models_reg_dev.c')

# prepare new model directory + files
def ptatinModelPrepare():
  dest = os.path.join('sandbox',mname)
  shell.copytree('template',dest)
  print('[MGenerator]: Creating dir ' + dest)

  src = os.path.join('sandbox',mname,'model_ops_template.c')
  dest = os.path.join('sandbox',mname,mname+'_model_def.c')
  shell.move(src,dest)

  #src = os.path.join('sandbox',mname,'model_template_ctx.h')
  #dest = os.path.join('sandbox',mname,mname+'_model_context.h')
  #shell.move(src,dest)


#ptatinModelPrepare()

ptatinModelRegUpdate()

