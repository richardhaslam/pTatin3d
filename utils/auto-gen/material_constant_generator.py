
import os
import socket
import sys
import time
import datetime



def write_out_c_class_externdefs( ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);
	if L != len(variable_type_list):
		print 'ERROR: length variable_name_list[] != variable_name_list[]'
	if L != len(variable_textural_name_list):
		print 'ERROR: length variable_textural_name_list[] != variable_name_list[]'
	

	#print '\n#error(<<REMOVE AUTOGENERATED TAG>> ================== FILE ['+ParticleClass+'_def.c] ==================)'

	# comment
	now = datetime.datetime.now()
	print '/*'
	print '  Auto generated by version 0.0 of material_constant_generator.py'
	print '  on '+ socket.gethostname() +', at '+str(now)+' by '+os.getenv('USER')
	print '*/\n'


	print '#include <stdio.h>'
	print '#include <string.h>'
	print '#include <stdlib.h>'
	print '#include <petsc.h>\n'
	print '#include "'+ParticleClass+'_def.h"\n\n'

	# function prototype
	print 'PetscErrorCode MaterialConstantsReportParseError(const char model_name[],const char field_name[],const int region);\n\n'

	# generate classname
	print 'const char '+ParticleClass+'_classname[] = "'+ParticleClass+'";\n'

	# strings and constant stuff #
	print 'const int ' + ParticleClass +'_nmembers = ' + str(L) + ';\n'
	
	
	# sizes #
	l = 'const size_t ' + ParticleClass +'_member_sizes[] = {'
	print l
	for f in xrange(L):
		l = '  ' + str(variable_extend_list[f])  + ' * sizeof(' + str(variable_type_list[f]) +')'
		if f < L-1:
			l = l + ','
		print l
	print '};\n'
	

	# string names #
	l = 'const char *' + ParticleClass +'_member_names[] = {'
	print l
	for f in xrange(L):
		l = '  ' + '"' + variable_textural_name_list[f] + '"'
		if f < L-1:
			l = l + ','
		print l
		
	l = '};\n'
	print l



def write_out_getters( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	# comment
	if prototype == 'False':
		print '\n/* ================================================================= */'
		print '/*   Getters for '+ParticleClass+' */'
		print '/* ================================================================= */'

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if prototype == 'True':
				print 'void '+ ParticleClass +'GetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' *data);'
				continue
			

			l = 'void '+ ParticleClass +'GetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' *data) \n{'
			print l

			l = '  *data = point->'+variable_name_list[f]+';'
			print l
			print '}\n'
		
		else:

			if prototype == 'True':
				print 'void '+ ParticleClass +'GetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' *data[]);'
				continue

			l = 'void '+ ParticleClass +'GetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' *data[]) \n{'
			print l


			l = '  *data = point->'+variable_name_list[f]+';'
			print l
			print '}\n'
		

def write_out_setters( protoype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if protoype == 'False':
		print '\n/* ================================================================= */'
		print '/*   Setters for '+ParticleClass+' */'
		print '/* ================================================================= */'

	for f in xrange(L):

		if variable_extend_list[f] == 1:
			if protoype == 'True':
				print 'void '+ ParticleClass +'SetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' data);'
				continue

			l = 'void '+ ParticleClass +'SetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' data) \n{'
			print l

			l = '  point->'+variable_name_list[f]+' = data;'
			print l
			print '}\n'
		
		else:
			if protoype == 'True':
				print 'void '+ ParticleClass +'SetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' data[]);'
				continue

			l = 'void '+ ParticleClass +'SetField_'+variable_textural_name_list[f]+'('+ParticleClass+' *point,'+variable_type_list[f]+' data[]) \n{'
			print l


			l = '  memcpy( &point->'+variable_name_list[f]+'[0], data, sizeof('+variable_type_list[f]+')*'+str(variable_extend_list[f])+' );'
			print l
			
			print '}\n'



# -------------------------------------------------------------------
def write_out_SetDefault( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);


	if protoype == 'True':
		l = 'PetscErrorCode MaterialConstantsSetDefault_'+ ClassNameShort +'(int nr,'+ ClassName +' _data[],'
		for f in xrange(L-1):
			l = l + variable_type_list[f] + ' ' +  variable_name_list[f] + ','
		l = l + variable_type_list[L-1] + ' ' +  variable_name_list[L-1] + ' );\n'
		print(l)
		return
	
	print '#undef __FUNCT__'
	print '#define __FUNCT__ \"MaterialConstantsSetDefault_'+ ClassNameShort+'\"'
	print 'PetscErrorCode MaterialConstantsSetDefault_'+ ClassNameShort +'( \n    int nr,'+ ClassName +' _data[],'
	for f in xrange(L-1):
		print '    ' + variable_type_list[f] + ' ' +  variable_name_list[f] + ','
	print '    ' + variable_type_list[L-1] + ' ' +  variable_name_list[L-1] + ' )\n{'


	print '  int r; \n'

	print '  for (r=0; r<nr; r++) {'
	
	for f in xrange(L):
		print '    _data[r].' + variable_name_list[f] + ' =  ' + variable_name_list[f] + ';'
	
	print '  }\n'

	print '  PetscFunctionReturn(0);'
	print '} \n'

def write_out_SetDefaultAll( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);


	if protoype == 'True':
		print 'void MaterialConstantsSetDefaultAll_'+ ClassNameShort +'(int nr,'+ ClassName +' _data[]);'
		return
	
	print 'void MaterialConstantsSetDefaultAll_'+ ClassNameShort +'( \n    int nr,'+ ClassName +' _data[])\n{'

	print '  int r; \n'

	print '  for (r=0; r<nr; r++) {'
	
#	for f in xrange(L):
#		print '    _data[r].' + variable_name_list[f] + ' =  ' + variable_name_list[f] + ';'

	for f in xrange(L):
		if variable_type_list[f] == 'float':
			print '    { double float;'
			print '      ' + ClassName + 'GetDefault_' + variable_textural_name_list[f] + '((float*)&value);'
			if variable_extend_list[f] == 1:
				print '      ' + ClassName + 'SetField_' + variable_textural_name_list[f] + '(&_data[r],(float)value);'
			else:
				for ii in xrange(variable_extend_list[f]):
					print '      _data[r].' + variable_name_list[f] + '[' + str(ii) + ']' + ' = value;'

			print '    }\n'

		if variable_type_list[f] == 'double':
			print '    { double value;'
			print '      ' + ClassName + 'GetDefault_' + variable_textural_name_list[f] + '((double*)&value);'
			if variable_extend_list[f] == 1:
				print '      ' + ClassName + 'SetField_' + variable_textural_name_list[f] + '(&_data[r],(double)value);'
			else:
				for ii in xrange(variable_extend_list[f]):
					print '      _data[r].' + variable_name_list[f] + '[' + str(ii) + ']' + ' = value;'

			print '    }\n'

		if variable_type_list[f] == 'int':
			print '    { int value;'
			print '      ' + ClassName + 'GetDefault_' + variable_textural_name_list[f] + '((int*)&value);'
			if variable_extend_list[f] == 1:
				print '      ' + ClassName + 'SetField_' + variable_textural_name_list[f] + '(&_data[r],(int)value);'
			else:
				for ii in xrange(variable_extend_list[f]):
					print '      _data[r].' + variable_name_list[f] + '[' + str(ii) + ']' + ' = value;'

			print '    }\n'

	
	print '  }\n'

	print '} \n'

def write_out_GetDefault( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults ):

	#check length
	L = len(variable_name_list);

	# comment
	if prototype == 'False':
		print '\n/* ================================================================= */'
		print '/*   Getters for default parameters ('+ParticleClass+') */'
		print '/* ================================================================= */'

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if prototype == 'True':
				print 'void '+ ParticleClass +'GetDefault_'+variable_textural_name_list[f]+'('+variable_type_list[f]+' *data);'
				continue
			

			l = 'void '+ ParticleClass +'GetDefault_'+variable_textural_name_list[f]+'('+variable_type_list[f]+' *data) \n{'
			print l

			l = '  *data = ' + '(' + variable_type_list[f] + ')' +str(variable_defaults[f])+';'
			print l
			print '}\n'
		
		else:

			if prototype == 'True':
				print 'void '+ ParticleClass +'GetDefault_'+variable_textural_name_list[f]+'('+variable_type_list[f]+' *data);'
				continue

			l = 'void '+ ParticleClass +'GetDefault_'+variable_textural_name_list[f]+'('+variable_type_list[f]+' *data) \n{'
			print l


			l = '  *data = '+ '(' + variable_type_list[f] + ')' +str(variable_defaults[f])+';'
			print l
			print '}\n'


# -------------------------------------------------------------------
def write_out_SetValues( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if protoype == 'True':
		l = 'PetscErrorCode MaterialConstantsSetValues_'+ ClassNameShort +'(const int region_id,'+ ClassName +' _data[],'
		for f in xrange(L-1):
			if variable_extend_list[f] == 1:
				l = l + variable_type_list[f] + ' ' +  variable_name_list[f] + ','
			else:
				l = l + variable_type_list[f] + ' *' +  variable_name_list[f] + ','

		if variable_extend_list[L-1] == 1:
			l = l + variable_type_list[L-1] + ' ' +  variable_name_list[L-1] + ');'
		else:
			l = l + variable_type_list[L-1] + ' *' +  variable_name_list[L-1] + ');'

		print(l)
		return


	print '#undef __FUNCT__'
	print '#define __FUNCT__ \"MaterialConstantsSetValues_'+ ClassNameShort+'\"'
	print 'PetscErrorCode MaterialConstantsSetValues_'+ ClassNameShort +'(const int region_id,'+ ClassName +' _data[],'
	for f in xrange(L-1):
		if variable_extend_list[f] == 1:
			print '    ' + variable_type_list[f] + ' ' +  variable_name_list[f] + ','
		else:
			print '    ' + variable_type_list[f] + ' *' +  variable_name_list[f] + ','

	if variable_extend_list[L-1] == 1:
		print '    ' + variable_type_list[L-1] + ' ' +  variable_name_list[L-1] + ')\n{'
	else:
		print '    ' + variable_type_list[L-1] + ' *' +  variable_name_list[L-1] + ')\n{'

	print '  ' + ClassName + ' *data = &_data[region_id];'

	for f in xrange(L):
		if variable_extend_list[f] == 1:
			print '  data->' + variable_name_list[f] + ' =  ' + variable_name_list[f] + ';'
		else:	
			print '  memcpy(&data->' + variable_name_list[f] + ',' + variable_name_list[f] + ',' + str(variable_extend_list[f]) + '*sizeof(' + variable_type_list[f] +')' +');'


	print '  PetscFunctionReturn(0);'
	print '} \n'

def write_out_ScaleValues( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length

	floatlist = []
	floatvarlist = []
	for item in xrange(len(variable_type_list)):
		if variable_type_list[item] != 'int':
			floatlist.append(variable_type_list[item])
			floatvarlist.append(variable_name_list[item])
	L = len(floatlist);
	

	if protoype == 'True':
		l = 'PetscErrorCode MaterialConstantsScaleValues_'+ ClassNameShort +'(const int region_id,'+ ClassName +' _data[],'
		for f in xrange(L-1):
			l = l + floatlist[f] + ' ' +  floatvarlist[f] + ','
		l = l + floatlist[L-1] + ' ' +  floatvarlist[L-1] + ');\n'

		print(l)
		return

	print '#undef __FUNCT__'
	print '#define __FUNCT__ \"MaterialConstantsScaleValues_'+ ClassNameShort+'\"'
	print 'PetscErrorCode MaterialConstantsScaleValues_'+ ClassNameShort +'(const int region_id,'+ ClassName +' _data[],'
	for f in xrange(L-1):
		print '    ' + floatlist[f] + ' ' +  floatvarlist[f] + ','
	print '    ' + floatlist[L-1] + ' ' +  floatvarlist[L-1] + ')\n{'


	print '  ' + ClassName + ' *data = &_data[region_id];\n'

	L = len(variable_name_list);
	for f in xrange(L):

		if variable_type_list[f] == 'float':
			print '  { double float;'
			print '    ' + ClassName + 'GetField_' + variable_textural_name_list[f] + '(data,(float*)&value);'
			print '    value = value / ' +  variable_name_list[f] + ';'
			print '    ' + ClassName + 'SetField_' + variable_textural_name_list[f] + '(data,(float)value);'
			print '  }\n'

		if variable_type_list[f] == 'double':
			print '  { double value;'
			print '    ' + ClassName + 'GetField_' + variable_textural_name_list[f] + '(data,(double*)&value);'
			print '    value = value / ' +  variable_name_list[f] + ';'
			print '    ' + ClassName + 'SetField_' + variable_textural_name_list[f] + '(data,(double)value);'
			print '  }\n'

	

	print '  PetscFunctionReturn(0);'
	print '} \n'

def write_out_SetFromOptions( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);


	if protoype == 'True':
		print 'PetscErrorCode MaterialConstantsSetFromOptions_'+ ClassNameShort +'(const char model_name[],const int region_id,' + ClassName +' _data[],PetscBool essential);'

		return

	type   = [ 'float', 'double', 'int' ]

	print '#undef __FUNCT__'
	print '#define __FUNCT__ \"MaterialConstantsSetFromOptions_'+ ClassNameShort+'\"'
	print 'PetscErrorCode MaterialConstantsSetFromOptions_'+ ClassNameShort +'(const char model_name[],const int region_id,' + ClassName +' _data[],PetscBool essential)\n{'

	print '  char                         opt_name[PETSC_MAX_PATH_LEN];'
	print '  PetscBool                    found;'
	print '  PetscErrorCode               ierr;\n'
	print '  ' + ClassName + ' *data = &_data[region_id];'

	for f in xrange(L):
		print '  /* options for ' + variable_textural_name_list[f] + ' ==>> ' + variable_name_list[f] + ' */'
		print '  sprintf(opt_name,"-' + variable_textural_name_list[f] + '_%d",region_id);'

		if variable_type_list[f] == 'float':
			if variable_extend_list[f] == 1:
				print '  { PetscReal value;'
				print '    ierr = PetscOptionsGetReal(NULL,model_name,opt_name,&value,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				print '      data->' + variable_name_list[f] + ' = (float)value;'
			else:
				print '  { PetscReal value['+ str(variable_extend_list[f]) + '];'
				print '    PetscInt nv;'
				print '    ierr = PetscOptionsGetRealArray(NULL,model_name,opt_name,value,&nv,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				for ii in xrange(variable_extend_list[f]):
					print '      data->' + variable_name_list[f] + '[' + str(ii) + ']' + ' = (float)value[' + str(ii) + '];'
			print '    }'

			print '    else if ( (!found)  && (essential) ) {'
			print '      ierr = MaterialConstantsReportParseError(model_name,"' + variable_textural_name_list[f] + '",region_id);CHKERRQ(ierr);'
			print '  }}\n'

		if variable_type_list[f] == 'double':
			if variable_extend_list[f] == 1:
				print '  { PetscReal value;'
				print '    ierr = PetscOptionsGetReal(NULL,model_name,opt_name,&value,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				print '      data->' + variable_name_list[f] + ' = (double)value;'
			else:
				print '  { PetscReal value['+ str(variable_extend_list[f]) + '];'
				print '    PetscInt nv;'
				print '    ierr = PetscOptionsGetRealArray(NULL,model_name,opt_name,value,&nv,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				for ii in xrange(variable_extend_list[f]):
					print '      data->' + variable_name_list[f] + '[' + str(ii) + ']' + ' = (double)value[' + str(ii) + '];'
			print '    }'

			print '    else if ( (!found)  && (essential) ) {'
			print '      ierr = MaterialConstantsReportParseError(model_name,"' + variable_textural_name_list[f] + '",region_id);CHKERRQ(ierr);'
			print '  }}\n'

		if variable_type_list[f] == 'int':
			if variable_extend_list[f] == 1:
				print '  { PetscInt value;'
				print '    ierr = PetscOptionsGetInt(NULL,model_name,opt_name,&value,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				print '      data->' + variable_name_list[f] + ' = (int)value;'
			else:
				print '  { PetscInt value['+ str(variable_extend_list[f]) + '];'
				print '    PetscInt nv;'
				print '    ierr = PetscOptionsGetIntArray(NULL,model_name,opt_name,value,&nv,&found);CHKERRQ(ierr);'
				print '    if (found) {'
				for ii in xrange(variable_extend_list[f]):
					print '      data->' + variable_name_list[f] + '[' + str(ii) + ']' + ' = (int)value[' + str(ii) + '];'
			print '    }'

			print '    else if ( (!found)  && (essential) ) {'
			print '      ierr = MaterialConstantsReportParseError(model_name,"' + variable_textural_name_list[f] + '",region_id);CHKERRQ(ierr);'
			print '  }}\n'


	print '  PetscFunctionReturn(0);'
	print '} \n'



def write_out_PrintValues( protoype, ClassName, ClassNameShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if protoype == 'True':
		print 'PetscErrorCode MaterialConstantsPrintValues_'+ ClassNameShort +'(const char model_name[],const int region_id,' + ClassName +' _data[]);'

		return


	print '#undef __FUNCT__'
	print '#define __FUNCT__ \"MaterialConstantsPrintValues_'+ ClassNameShort+'\"'
	print 'PetscErrorCode MaterialConstantsPrintValues_'+ ClassNameShort +'(const char model_name[],const int region_id,' + ClassName +' _data[]) \n{'
	print '  ' + ClassName + ' *data = &_data[region_id];'
	print '  char   opt_name[PETSC_MAX_PATH_LEN];\n'


	print '  PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------------------------------------------------------\\n");'
	print '  PetscPrintf(PETSC_COMM_WORLD,"  MaterialView(' + ClassNameShort + '): RegionIndex[%d]\\n", region_id);'
	for f in xrange(L):
		print '  /* options for ' + variable_textural_name_list[f] + ' ==>> ' + variable_name_list[f] + ' */'
		print '  sprintf(opt_name,"-%s_'+variable_textural_name_list[f]+'_%d", model_name,region_id);'

		if variable_extend_list[f] == 1:

			if variable_type_list[f] == 'float':
				print '  { float value;'
				print '    ' + ClassName + 'GetField_' + variable_textural_name_list[f]+ '(data,(float*)&value);'
				print '    PetscPrintf(PETSC_COMM_WORLD,"    ' + variable_textural_name_list[f] + ' = %1.4e (%s) \\n", value,opt_name); \n  }\n'

			if variable_type_list[f] == 'double':
				print '  { double value;'
				print '    ' + ClassName + 'GetField_' + variable_textural_name_list[f]+ '(data,(double*)&value);'
				print '    PetscPrintf(PETSC_COMM_WORLD,"    ' + variable_textural_name_list[f] + ' = %1.4e (%s) \\n", value,opt_name); \n  }\n'

			if variable_type_list[f] == 'int':
				print '  { int value;'
				print '    ' + ClassName + 'GetField_' + variable_textural_name_list[f]+ '(data,(int*)&value);'
				print '    PetscPrintf(PETSC_COMM_WORLD,"    ' + variable_textural_name_list[f] + ' = %d (%s) \\n", value,opt_name); \n  }\n'

	
	print '  PetscFunctionReturn(0);'
	print '} \n'


def write_out_viewer( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if prototype == 'False':
		print '\n/* ================================================================= */'
		print '/*   C-viewer for '+ParticleClass+' */'
		print '/* ================================================================= */'

	if prototype == 'True':
		print 'void '+ParticleClass+'View('+ParticleClass+' *point);'
		return

	l = 'void '+ParticleClass+'View('+ParticleClass+' *point)\n{'
	print l
	
	type   = [ 'float', 'double', 'char', 'int', 'long int', 'short' ]
	format = [ '%1.6e', '%1.6e' ,   '%c',  '%d', '%ld',      '%d' ]
	
	
	for f in xrange(L):

		print '  {'
		if variable_extend_list[f] == 1:
			print '    '+variable_type_list[f]+' data;'
		else:
			print '    '+variable_type_list[f]+' *data;'
		
		print '    '+ParticleClass+'GetField_'+variable_textural_name_list[f]+'(point,&data);' 


		if variable_extend_list[f] == 1:

			format_i = ''
			for i in xrange(len(type)):
				if variable_type_list[f] == type[i]:
					format_i = format[i]
					break
			
			print '    printf(\"field: '+variable_textural_name_list[f]+' = '+format_i+'; [size %zu; type '+variable_type_list[f]+'; variable_name '+variable_name_list[f]+']\\n\",data, '+ParticleClass+'_member_sizes['+str(f)+'] );'
			
		else:

			format_i = ''
			for i in xrange(len(type)):
				if variable_type_list[f] == type[i]:
					format_i = format[i]
					break

			for K in xrange(variable_extend_list[f]):
				print '    printf(\"field: '+variable_textural_name_list[f]+'['+str(K)+'] = '+format_i+'; [size %zu; type '+variable_type_list[f]+'; variable_name '+variable_name_list[f]+']\\n\",data['+str(K)+'], '+ParticleClass+'_member_sizes['+str(f)+'] );'


		# close stack fram
		print '  }'


	print '}\n'
	
def write_out_c_class( ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults ):

	#check length
	L = len(variable_name_list);
	if L != len(variable_type_list):
		print 'ERROR: length variable_name_list[] != variable_name_list[]'
	if L != len(variable_textural_name_list):
		print 'ERROR: length variable_textural_name_list[] != variable_name_list[]'
	

	#print '\n#error(<<REMOVE AUTOGENERATED TAG>> ================== FILE ['+ParticleClass+'_def.c] ==================)'

	# comment
	now = datetime.datetime.now()
	print '/*'
	print '  Auto generated by version 0.0 of material_constant_generator.py'
	print '  on '+ socket.gethostname() +', at '+str(now)+' by '+os.getenv('USER')
	print '*/\n'

	
	# safe macro
	print '#ifndef __'+ParticleClass+'_DEF_H__'
	print '#define __'+ParticleClass+'_DEF_H__\n'


	# generate struct
	print 'typedef struct {'
	for f in xrange(L):
		if variable_extend_list[f] == 1:
			print(' ' + variable_type_list[f] + ' ' + variable_name_list[f] + ';')
		else:
			print(' ' + variable_type_list[f] + ' ' + variable_name_list[f] + '[' + str(variable_extend_list[f]) + '];')

	print('} ' + ParticleClass + ';\n')

	# enums #
	print 'typedef enum {'
	for f in xrange(L):
		if f == 0:
			l = '  ' + ParticleClassShortName + '_' + variable_textural_name_list[f] + ' = 0'
		else:
			l = '  ' + ParticleClassShortName + '_' + variable_textural_name_list[f]

		if f < L-1:
			l = l + ','
		print l
		
	print('} ' + ParticleClass + 'TypeName;\n')

	print 'extern const char '+ParticleClass+'_classname[];\n'

	print 'extern const int ' + ParticleClass +'_nmembers;\n'
	
	print 'extern const size_t ' + ParticleClass +'_member_sizes[];\n'

	print 'extern const char *' + ParticleClass +'_member_names[];\n'


	# dump prototypes here
	print '/* prototypes */'
	write_out_getters( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_setters( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_viewer(  'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


	print ''
	#write_out_SetDefault( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_GetDefault( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults )
	write_out_SetDefaultAll( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_out_SetFromOptions( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_PrintValues( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_SetValues('True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_ScaleValues('True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	print '\n#endif'


def MATERIALPROP_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults ):

	# write out the header
	file = open(ParticleClass+'_def.h','w')
	sys.stdout = file
	
	write_out_c_class( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults )

	file.close()
	sys.stdout = sys.__stdout__



	# write out the c file
	file = open(ParticleClass+'_def.c','w')
	sys.stdout = file

	write_out_c_class_externdefs( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_out_getters( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_setters( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_viewer( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

#	write_out_SetDefault( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_GetDefault( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list, variable_defaults )
	write_out_SetDefaultAll( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_out_SetFromOptions( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_PrintValues( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_SetValues('False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_ScaleValues('False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	file.close()
	sys.stdout = sys.__stdout__



# ==================================================================================
# Material constants 


def Generate_pTatin_MaterialConst_ViscosityArrh():
	ClassName      = 'MaterialConst_ViscosityArrh'
	ClassNameShort = 'ViscosityArrh'
	variable_names =          [ 'preexpA','Ascale' ,'entalpy' , 'Vmol'   ,'nexp'    ,'Tref'    ,'Eta_scale','P_scale']       
	variable_types =          [ 'double' ,'double' , 'double' , 'double' , 'double' , 'double' ,'double'   , 'double' ]
	variable_extents        = [ 1        ,      1  ,       1  , 1        , 1        ,  1       , 1         ,    1    ]
	variable_textural_names = [ 'preexpA','Ascale' ,'entalpy' , 'Vmol'   ,'nexp'    ,'Tref'    ,'Eta_scale','P_scale']       
	variable_defaults       = [ 1.0      , 1.0     , 0.0      , 0.0      , 0        ,  0.0     , 1.0       , 1.0     ]

	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_pTatin_MaterialConst_DiffusivityConst():
	ClassName      = 'MaterialConst_DiffusivityConst'
	ClassNameShort = 'DiffusivityConst'
	variable_names =          [ 'k',      'Cp',     'rho_t'  ]       
	variable_types =          [ 'double', 'double', 'double' ]
	variable_extents        = [ 1,         1,       1        ]
	variable_textural_names = [ 'conductivity','heat_capacity' ,'thermal_density' ]       
	variable_defaults       = [ 1.0,       1.0,      1.0     ]

	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )

def Generate_pTatin_MaterialConst_VolumetricHeatingConst():
	ClassName      = 'MaterialConst_VolumetricHeatingConst'
	ClassNameShort = 'VolumetricHeatingConst'
	variable_names =          [ 'H',     ]       
	variable_types =          [ 'double' ]
	variable_extents        = [ 1        ]
	variable_textural_names = [ 'heat_source' ]       
	variable_defaults       = [ 0.0   ]

	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )

def Generate_pTatin_MaterialConst_PowerLaw():
	ClassName      = 'MaterialConst_PowerLaw'
	ClassNameShort = 'PowerLaw'
	variable_names =          [ 'A'        , 'n'      ]       
	variable_types =          [ 'double'   , 'double' ]
	variable_extents        = [ 1          , 1        ]
	variable_textural_names = [ 'prefactor', 'n'      ]
	variable_defaults       = [ 0.0        , 0.0      ]

	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


# Material constants for energy equation
def Generate_EnergyMaterialConstants():
	ClassName      = 'EnergyMaterialConstants'
	ClassNameShort = 'EnergyMaterialConstants'
	variable_names =          [ 'alpha', 'beta', 'rho_ref', 'Cp', 'density_type', 'conductivity_type', 'source_type'   ]
	variable_types =          [ 'double', 'double', 'double', 'double', 'int', 'int', 'int' ]
	variable_extents        = [  1,       1,         1,        1,        1,1,7    ]
	variable_textural_names = [ 'ThermalExpansivity', 'Compressibility', 'ReferenceDensity', 'SpecificHeat', 'DensityMethod', 'ConductivityMethod', 'SourceMethod' ]
	variable_defaults       = [ 0.0, 0.0, 0.0, 0.0, 0,0,0 ]
	
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_EnergySourceConst():
	ClassName      = 'EnergySourceConst'
	ClassNameShort = 'SourceConst'
	variable_names =          [ 'H'          ]
	variable_types =          [ 'double'     ]
	variable_extents        = [ 1            ]
	variable_textural_names = [ 'HeatSource' ]
	variable_defaults       = [ 0.0          ]
        
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_EnergySourceDecay():
	ClassName      = 'EnergySourceDecay'
	ClassNameShort = 'SourceDecay'
	variable_names =          [ 'H0', 'lambda'      ]
	variable_types =          [ 'double' , 'double' ]
	variable_extents        = [ 1, 1           ]
	variable_textural_names = [ 'HeatSourceRef', 'HalfLife' ]
	variable_defaults       = [ 0.0, 0.0 ]
	
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_EnergySourceAdiabaticAdvection():
	ClassName      = 'EnergySourceAdiabaticAdvection'
	ClassNameShort = 'SourceAdiabaticAdv'
	variable_names =          [ 'dTdy'   ]
	variable_types =          [ 'double' ]
	variable_extents        = [ 1        ]
	variable_textural_names = [ 'VerticalThermalGradient' ]
	variable_defaults       = [ 0.0 ]
	
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_EnergyConductivityConst():
	ClassName      = 'EnergyConductivityConst'
	ClassNameShort = 'ConductivityConst'
	variable_names =          [ 'k0' ]
	variable_types =          [ 'double' ]
	variable_extents        = [ 1 ]
	variable_textural_names = [ 'k0' ]
	variable_defaults       = [ 0.0 ]
	
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )


def Generate_EnergyConductivityThreshold():
	ClassName      = 'EnergyConductivityThreshold'
	ClassNameShort = 'ConductivityThreshold'
	variable_names =          [ 'k0', 'k1', 'T_threshold', 'dT'     ]
	variable_types =          [ 'double' , 'double', 'double', 'double' ]
	variable_extents        = [ 1, 1, 1, 1           ]
	variable_textural_names = [ 'k0', 'k1', 'ThresholdTemperature', 'DeltaT' ]
	variable_defaults       = [ 0.0, 0.0, 0.0, 0.0 ]
	
	MATERIALPROP_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names, variable_defaults )





# Call all functions to generate all data types
def main():
	
	## material constants ##
	Generate_pTatin_MaterialConst_ViscosityArrh()
	Generate_pTatin_MaterialConst_DiffusivityConst()
	Generate_pTatin_MaterialConst_VolumetricHeatingConst()
	Generate_pTatin_MaterialConst_PowerLaw()

	## Material metadata for energy equations
	Generate_EnergyMaterialConstants()
	Generate_EnergySourceConst()
	Generate_EnergySourceDecay()
	Generate_EnergySourceAdiabaticAdvection()

	Generate_EnergyConductivityConst()
	Generate_EnergyConductivityThreshold()




if __name__ == "__main__":
	main()


