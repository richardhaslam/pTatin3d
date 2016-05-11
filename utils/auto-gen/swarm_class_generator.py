
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
	

	print '\n#error(<<REMOVE AUTOGENERATED TAG>> ================== FILE ['+ParticleClass+'_def.c] ==================)'

	# comment
	now = datetime.datetime.now()
	print '/*'
	print '  Auto generated by version 0.0 of swarm_class_generator.py'
	print '  on '+ socket.gethostname() +', at '+str(now)+' by '+os.getenv('USER')
	print '*/\n'


	print '#include <stdio.h>'
	print '#include <string.h>'
	print '#include <mpi.h>\n'
	print '#include "'+ParticleClass+'_def.h"\n\n'

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

	# MPI type #
	print('MPI_Datatype MPI_' + ParticleClass.upper() + ';\n')
	


def write_out_getters( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	# comment
	if prototype == 'False':
		print '\n/* ===================================== */'
		print '/* Getters for '+ParticleClass+' */'
		print '/* ===================================== */'

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
		print '\n/* ===================================== */'
		print '/* Setters for '+ParticleClass+' */'
		print '/* ===================================== */'

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



def write_out_viewer( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if prototype == 'False':
		print '\n/* ===================================== */'
		print '/* C-viewer for '+ParticleClass+' */'
		print '/* ===================================== */'

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
	


def write_vtk_viewer( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if prototype == 'False':
		print '\n/* ===================================== */'
		print '/* VTK viewer for '+ParticleClass+' */'
		print '/* ===================================== */'
	else:
		print 'void '+ ParticleClass +'VTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const '+ParticleClass+' points[]);'
		return
	

	l = 'void '+ ParticleClass +'VTKWriteAsciiAllFields(FILE *vtk_fp,const int N,const '+ParticleClass+' points[]) \n{'
	print l
	print '  int p;'

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if variable_type_list[f] == 'double':
				vtk_type = 'Float64'
				cast     = '(double)'
				format   = '%lf'
			elif variable_type_list[f] == 'float':
				vtk_type = 'Float32'
				cast     = '(float)'
				format   = '%f'
			elif variable_type_list[f] == 'int':
				vtk_type = 'Int32'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'long int':
				vtk_type = 'Int64'
				cast     = '(long int)'
				format   = '%ld'
			elif variable_type_list[f] == 'char':
				vtk_type = 'Int8'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'short':
				vtk_type = 'Int16'
				cast     = '(short)'
				format   = '%d'
			else:
				print 'Unknown type: Cannot find equivalent VTKType'
				exit()

			##
			## maybe here you want to use variable_textural_name_list rather than variable_name_list[] ??
			##
			l = '  fprintf( vtk_fp, "\\t\\t\\t\\t<DataArray type=\\"'+vtk_type+'\\" Name=\\"'+variable_name_list[f]+'\\" format=\\"ascii\\">\\n");'
			print l

			print '  for(p=0;p<N;p++) {'
			print '    fprintf( vtk_fp,"\\t\\t\\t\\t\\t'+format+'\\n",'+cast+'points[p].'+variable_name_list[f]+');'
			print '  }'

			print '  fprintf( vtk_fp, "\\t\\t\\t\\t</DataArray>\\n");'

	print '}\n'


def write_vtk_viewer_binary_appended_header( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if prototype == 'False':
		print '\n/* ===================================== */'
		print '/* VTK binary (appended header) viewer for '+ParticleClass+' */'
		print '/* ===================================== */'
	else:
		print 'void '+ ParticleClass +'VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const '+ParticleClass+' points[]);'
		return
	

	l = 'void '+ ParticleClass +'VTKWriteBinaryAppendedHeaderAllFields(FILE *vtk_fp,int *offset,const int N,const '+ParticleClass+' points[]) \n{'
	print l
	# declare variables here

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if variable_type_list[f] == 'double':
				vtk_type = 'Float64'
				cast     = '(double)'
				format   = '%lf'
			elif variable_type_list[f] == 'float':
				vtk_type = 'Float32'
				cast     = '(float)'
				format   = '%f'
			elif variable_type_list[f] == 'int':
				vtk_type = 'Int32'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'long int':
				vtk_type = 'Int64'
				cast     = '(long int)'
				format   = '%ld'
			elif variable_type_list[f] == 'char':
				vtk_type = 'Int8'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'short':
				vtk_type = 'Int16'
				cast     = '(short)'
				format   = '%d'
			else:
				print 'Unknown type: Cannot find equivalent VTKType'
				exit()

			##
			## maybe here you want to use variable_textural_name_list rather than variable_name_list[] ??
			##
			l = '  fprintf( vtk_fp, "\\t\\t\\t\\t<DataArray type=\\"'+vtk_type+'\\" Name=\\"'+variable_name_list[f]+'\\" format=\\"appended\\"  offset=\\"%d\\" />\\n",*offset);'
			print l
			l = '  *offset = *offset + sizeof(int) + N * sizeof(' + variable_type_list[f] + ');\n'
			print l
		
		else:
			print '  /* Warning: swarm_class_generator.py is ignoring multi-component field ' + variable_name_list[f] + '[] */\n'
	
	print '}\n'

def write_vtk_viewer_binary_appended_data( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);

	if prototype == 'False':
		print '\n/* ================================================== */'
		print '/* VTK binary (appended data) viewer for '+ParticleClass+' */'
		print '/* ==================================================== */'
	else:
		print 'void '+ ParticleClass +'VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const '+ParticleClass+' points[]);'
		return

	l = 'void '+ ParticleClass +'VTKWriteBinaryAppendedDataAllFields(FILE *vtk_fp,const int N,const '+ParticleClass+' points[]) \n{'
	print l

	# declare variables here
	print '  int p,length;'
	print '  size_t atomic_size;\n'

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if variable_type_list[f] == 'double':
				vtk_type = 'Float64'
				cast     = '(double)'
				format   = '%lf'
			elif variable_type_list[f] == 'float':
				vtk_type = 'Float32'
				cast     = '(float)'
				format   = '%f'
			elif variable_type_list[f] == 'int':
				vtk_type = 'Int32'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'long int':
				vtk_type = 'Int64'
				cast     = '(long int)'
				format   = '%ld'
			elif variable_type_list[f] == 'char':
				vtk_type = 'Int8'
				cast     = '(int)'
				format   = '%d'
			elif variable_type_list[f] == 'short':
				vtk_type = 'Int16'
				cast     = '(short)'
				format   = '%d'
			else:
				print 'Unknown type: Cannot find equivalent VTKType'
				exit()


			l = '  atomic_size = sizeof(' + variable_type_list[f] + ');'
			print l
			print '  length = (int)( atomic_size * ((size_t)N) );'

			# write out length into file
			print '  fwrite( &length,sizeof(int),1,vtk_fp);'

			print '  for(p=0;p<N;p++) {'
			print '    fwrite( &points[p].'+variable_name_list[f]+',atomic_size,1,vtk_fp);'
			print '  }\n'

		else:
			print '  /* Warning: swarm_class_generator.py is ignoring multi-component field ' + variable_name_list[f] + '[] */\n'

	# close function
	print '}\n'


# Not sure how to treat the coordinates yet...
# <PPoints> ... </PPoints>
# are missing at the moment
#
def write_pvtu_viewerPPointData( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);


	if prototype == 'False':
		print '\n/* ===================================== */'
		print '/* PVTU viewer for '+ParticleClass+' */'
		print '/* ===================================== */'
	else:
		print 'void '+ ParticleClass +'PVTUWriteAllPPointDataFields(FILE *vtk_fp);'
		return
	

	l = 'void '+ ParticleClass +'PVTUWriteAllPPointDataFields(FILE *vtk_fp) \n{'
	print l

	# open point data
	#print '  fprintf(vtk_fp, "\\t\\t<PPointData>\\n");'

	for f in xrange(L):

		if variable_extend_list[f] == 1:

			if variable_type_list[f] == 'double':
				vtk_type = 'Float64'
			elif variable_type_list[f] == 'float':
				vtk_type = 'Float32'
			elif variable_type_list[f] == 'int':
				vtk_type = 'Int32'
			elif variable_type_list[f] == 'long int':
				vtk_type = 'Int64'
			elif variable_type_list[f] == 'char':
				vtk_type = 'Int8'
			elif variable_type_list[f] == 'short':
				vtk_type = 'Int16'
				cast     = '(short)'
				format   = '%d'
			else:
				print 'Unknown type: Cannot find equivalent VTKType'
				exit()


			print '  fprintf(vtk_fp, "\\t\\t\\t<PDataArray type=\\"'+ vtk_type +'\\" Name=\\"'+ variable_name_list[f] +'\\" NumberOfComponents=\\"'+ str(variable_extend_list[f]) +'\\"/>\\n");'

	# close point data
	#print '  fprintf(vtk_fp, "\\t\\t</PPointData>\\n");'


	print '}\n'


def write_out_mpi_type( prototype, ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):
    
	#check length
	L = len(variable_name_list);
    
	if prototype == 'False':
		print('\n/* ===================================== */')
		print('/* MPI data type for '+ParticleClass+' */')
		print('/* ===================================== */')
    
	if prototype == 'True':
		print 'int '+ParticleClass+'CreateMPIDataType(MPI_Datatype *ptype);'
		return
    
	l = 'int '+ParticleClass+'CreateMPIDataType(MPI_Datatype *ptype)\n{'
	print(l)

	l = '  MPI_Datatype newtype;'
	print(l)

	l = '  MPI_Datatype types[] = { '
	for f in xrange(L):
		if variable_type_list[f] == 'double':
			mpi_type = 'MPI_DOUBLE'
		elif variable_type_list[f] == 'float':
			mpi_type = 'MPI_FLOAT'
		elif variable_type_list[f] == 'int':
			mpi_type = 'MPI_INT'
		elif variable_type_list[f] == 'long int':
			mpi_type = 'MPI_LONG'
		elif variable_type_list[f] == 'char':
			mpi_type = 'MPI_CHAR'
		elif variable_type_list[f] == 'short':
			mpi_type = 'MPI_SHORT'
		else:
			print 'Unknown type: Cannot find equivalent MPI_Datatype'
			exit()
		
		l = l + mpi_type
		if f != (L-1):
			l = l + ' , '

	l = l + ' };'
	print(l)
	
	
	l = '  int blocklens[] = { '
	for f in xrange(L):
		l = l + str(variable_extend_list[f])
		if f != (L-1):
			l = l + ' , '
	l = l + ' };'
	print(l)


	l =     '  MPI_Aint loc[' + str(L+1) + '];\n'
	l = l + '  MPI_Aint disp[' + str(L) + '];\n'
	l = l + '  ' + ParticleClass + ' dummy;\n'
	l = l + '  int i,ierr;\n'
	print(l)

	l = '  ierr = MPI_Get_address(&dummy,&loc[0]);\n'
	for f in xrange(L):
		l = l + '  ierr = MPI_Get_address(&dummy.' + variable_name_list[f] + ',&loc[' + str(f+1) + ']);\n'
	print(l)
	

	l =     '  for (i=0; i<' + str(L) + '; i++) {\n'
	l = l + '    disp[i] = loc[i+1] - loc[0];\n'
	l = l + '  }\n'
	print(l)


	l =     '  ierr = MPI_Type_create_struct(' + str(L) + ',blocklens,disp,types,&newtype);\n'
	l = l + '  ierr = MPI_Type_commit(&newtype);'
	print(l)

	l =     '  *ptype = newtype;\n'
	l = l + '  return 0;\n'
	l = l + '}\n'
	print(l)


def write_out_c_class( ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	#check length
	L = len(variable_name_list);
	if L != len(variable_type_list):
		print 'ERROR: length variable_name_list[] != variable_name_list[]'
	if L != len(variable_textural_name_list):
		print 'ERROR: length variable_textural_name_list[] != variable_name_list[]'
	

	print '\n#error(<<REMOVE AUTOGENERATED TAG>> ================== FILE ['+ParticleClass+'_def.c] ==================)'

	# comment
	now = datetime.datetime.now()
	print '/*'
	print '  Auto generated by version 0.0 of swarm_class_generator.py'
	print '  on '+ socket.gethostname() +', at '+str(now)+' by '+os.getenv('USER')
	print '*/\n'

	
	# safe macro
	print '#ifndef __'+ParticleClass+'_DEF_H__'
	print '#define __'+ParticleClass+'_DEF_H__\n'


	print '#include <mpi.h>\n'

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
#	l = '} ' + ParticleClass + 'TypeName ;'
#	print l
#	print '\n'


	# generate classname
	#print 'const char '+ParticleClass+'_classname[] = "'+ParticleClass+'";\n'
	print 'extern const char '+ParticleClass+'_classname[];\n'

	

	# strings and constant stuff #
#	l = 'const int ' + ParticleClass +'_nmembers = ' + str(L) + ';'
	print 'extern const int ' + ParticleClass +'_nmembers;\n'
	
	
	# sizes #
#	l = 'const size_t ' + ParticleClass +'_member_sizes[] = {'
#	print l
#	for f in xrange(L):
#		l = '  ' + str(variable_extend_list[f])  + ' * sizeof(' + str(variable_type_list[f]) +')'
#		if f < L-1:
#			l = l + ','
#		print l
#	print '} ;'
#	print '\n'
	
	print 'extern const size_t ' + ParticleClass +'_member_sizes[];\n'

	# string names #
#	l = 'const char *' + ParticleClass +'_member_names[] = {'
#	print l
#	for f in xrange(L):
#		l = '  ' + '"' + variable_textural_name_list[f] + '"'
#		if f < L-1:
#			l = l + ','
#		print l
#		
#	l = '} ;'
#	print l
#	print '\n'

	print 'extern const char *' + ParticleClass +'_member_names[];\n'

	# MPI type #
	print('extern MPI_Datatype MPI_' + ParticleClass.upper() + ';\n')


	# dump prototypes here
	print '/* prototypes */'
	write_out_getters( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_setters( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_viewer( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_vtk_viewer( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_pvtu_viewerPPointData( 'True',ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


	write_vtk_viewer_binary_appended_header( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_vtk_viewer_binary_appended_data( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_out_mpi_type( 'True', ParticleClass, ParticleClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


	print '\n#endif'


def PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list ):

	# write out the header
	file = open(ParticleClass+'_def.h','w')
	sys.stdout = file
	
	write_out_c_class( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	file.close()
	sys.stdout = sys.__stdout__



	# write out the c file
	file = open(ParticleClass+'_def.c','w')
	sys.stdout = file

	write_out_c_class_externdefs( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_getters( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_setters( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_out_viewer( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_vtk_viewer( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_pvtu_viewerPPointData( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_vtk_viewer_binary_appended_header( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )
	write_vtk_viewer_binary_appended_data( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	write_out_mpi_type( 'False', ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

	file.close()
	sys.stdout = sys.__stdout__

# ==================================================================================
# Quadrature point fields

def Generate_pTatin3d_QuadraturePointVolumeCoefficientStokes():
	ParticleClass      = 'QPntVolCoefStokes'
	ParticleClassShort = 'QPVCStk'
	variable_name_list = [ 'eta',    'rho',    'Fu',     'Fp'     ]
	variable_type_list = [ 'double', 'double', 'double', 'double' ]
	variable_extend_list        = [ 1, 1, 3, 1 ]
	variable_textural_name_list = [ 'eta_effective', 'rho_effective', 'momentum_rhs', 'continuity_rhs' ]

	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

def Generate_pTatin3d_QuadraturePointSurfaceCoefficientStokes():
	ParticleClass      = 'QPntSurfCoefStokes'
	ParticleClassShort = 'QPSCStk'
	variable_name_list = [ 'normal', 'tangent1', 'tangent2', 'traction', 'eta','rho' ]
	variable_type_list = [ 'double', 'double',  'double', 'double', 'double', 'double'   ]
	variable_extend_list        = [ 3, 3, 3, 3, 1, 1 ]
	variable_textural_name_list = [ 'surface_normal', 'surface_tangent1', 'surface_tangent2', 'surface_traction', 'viscosity', 'density' ]


	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

def Generate_pTatin_QuadraturePointVolumeCoefficientEnergy():
	ParticleClass      = 'QPntVolCoefEnergy'
	ParticleClassShort = 'QPVCEgy'
	variable_name_list = [ 'diffusivity',    'heat_source' ]
	variable_type_list = [ 'double',         'double'      ]
	variable_extend_list        = [ 1, 1  ]
	variable_textural_name_list = [ 'diffusivity', 'heat_source' ]

	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


def Generate_pTatin_QuadraturePointVolumeCoefficientSPM_a():
	ClassName      = 'QPntVolCoefSPM'
	ClassShortName = 'QPVCSPM'
	variable_name_list = [ 'eff_diff_kappa',    'eff_source' ]
	variable_type_list = [ 'double',         'double'      ]
	variable_extend_list        = [ 1, 1  ]
	variable_textural_name_list = [ 'effective_diffusivity', 'effective_source' ]

	PARTICLE_CLASS_GENERATOR( ClassName, ClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

def Generate_pTatin_QuadraturePointVolumeCoefficientSPM_standard_fields():
	ClassName      = 'QPntVolCoefSPMStdFields'
	ClassShortName = 'QPVCSPMStdFields'
	variable_name_list = [ 'rho_r','rho_s', 'U_x1', 'U_x2', 'U' ]
	variable_type_list = [ 'double','double', 'double','double','double'  ]
	variable_extend_list        = [ 1,1, 1,1,1 ]
	variable_textural_name_list = [ 'density_rock', 'density_sediment', 'vel_horiz1','vel_horiz2','vel_vert' ]

	PARTICLE_CLASS_GENERATOR( ClassName, ClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

def Generate_pTatin_QuadraturePointVolumeCoefficientSPM_nonlinear_diffusivity():
	ClassName      = 'QPntVolCoefSPMNonlinearDiff'
	ClassShortName = 'QPVCSPMNonlinearDiff'
	variable_name_list = [ 'K','Sc', 'grad_z' ]
	variable_type_list = [ 'double','double', 'double'  ]
	variable_extend_list        = [ 1,1, 3 ]
	variable_textural_name_list = [ 'K', 'crticial_slope', 'slope' ]

	PARTICLE_CLASS_GENERATOR( ClassName, ClassShortName, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


# ==================================================================================
# Material point fields

def Generate_pTatin3d_MaterialPointStandard():
	ParticleClass      = 'MPntStd'
	ParticleClassShort = 'MPStd'
	variable_name_list = [ 'pid',      'coor',   'xi',     'phase', 'wil' ]
	variable_type_list = [ 'long int', 'double', 'double', 'int',   'int' ]
	variable_extend_list        = [ 1, 3, 3, 1, 1 ]
	variable_textural_name_list = [ 'point_index', 'global_coord', 'local_coord', 'phase_index', 'local_element_index' ]


	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )

def Generate_pTatin_MaterialPointStokesData():
	ParticleClass      = 'MPntPStokes'
	ParticleClassShort = 'MPPStk'
	variable_name_list = [ 'eta',      'rho'    ]
	variable_type_list = [ 'double',   'double' ]
	variable_extend_list        = [ 1, 1 ]
	variable_textural_name_list = [ 'eta_effective', 'density' ]


	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


def Generate_pTatin_MaterialPointStokesPlastic():
	ParticleClass      = 'MPntPStokesPl'
	ParticleClassShort = 'MPPStkPl'
	variable_name_list = [ 'e_plastic', 'is_yielding'    ]
	variable_type_list = [ 'float',   'short' ]
	variable_extend_list        = [ 1, 1 ]
	variable_textural_name_list = [ 'plastic_strain', 'yield_indicator' ]


	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


def Generate_pTatin_MaterialPointEnergy():
	ParticleClass      = 'MPntPEnergy'
	ParticleClassShort = 'MPPEgy'
	variable_name_list = [ 'diffusivity',    'heat_source' ]
	variable_type_list = [ 'double',         'double'      ]
	variable_extend_list        = [ 1, 1  ]
	variable_textural_name_list = [ 'diffusivity', 'heat_source' ]

	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


def Generate_pTatin_MaterialPointViscoElasticity():
	ParticleClass      = 'MPntPStokesVE'
	ParticleClassShort = 'MPPStkVE'
	variable_name_list = [ 'tau',      'mu' ]
	variable_type_list = [ 'double',   'double'      ]
	variable_extend_list        = [ 6, 1  ]
	variable_textural_name_list = [ 'deviatoric_stress', 'shear_modulus' ]

	PARTICLE_CLASS_GENERATOR( ParticleClass, ParticleClassShort, variable_name_list, variable_type_list, variable_extend_list, variable_textural_name_list )


# ==================================================================================
# Material constants 

def Generate_pTatin_MaterialConst_ViscosityConst():
	ClassName      = 'MaterialConst_ViscosityConst'
	ClassNameShort = 'ViscConst'
	variable_names =          [ 'eta0'   ]
	variable_types =          [ 'double' ]
	variable_extents        = [ 1        ]
	variable_textural_names = [ 'eta0'   ]
        
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )


def Generate_pTatin_MaterialConst_PlasticMises():
	ClassName      = 'MaterialConst_PlasticMises'
	ClassNameShort = 'PlasticMises'
	variable_names =          [ 'tau_yield'   , 'tau_yield_inf'    ]
	variable_types =          [ 'double'      , 'double'           ]
	variable_extents        = [ 1             , 1                  ]
	variable_textural_names = [ 'yield_stress', 'yield_stress_inf' ]
    
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_ViscosityZ():
	ClassName      = 'MaterialConst_ViscosityZ'
	ClassNameShort = 'ViscZ'
	variable_names =          [ 'eta0'      ,'zeta'     ,'zref'     ]
	variable_types =          [ 'double'    ,'double'   ,'double'   ]
	variable_extents        = [ 1           , 1         ,1          ]
	variable_textural_names = [ 'eta0'      ,'zeta'     ,'zref'     ]
    
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_PlasticDP():
	ClassName      = 'MaterialConst_PlasticDP'
	ClassNameShort = 'PlasticDP'
	variable_names =          [ 'phi'       ,'Co'       ,'phi_inf'      ,'Co_inf'       ,'tens_cutoff'  ,'hst_cutoff']
	variable_types =          [ 'double'    ,'double'   ,'double'       ,'double'       ,'double'       ,'double'    ]
	variable_extents        = [ 1           , 1         ,1              ,1              , 1             ,   1        ]
	variable_textural_names = [ 'friction'  ,'cohesion' ,'friction_inf' ,'cohesion_inf' ,'tens_cutoff'  ,'hst_cutoff']
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_MaterialType():
	ClassName      = 'MaterialConst_MaterialType'
	ClassNameShort = 'MaterialType'
	variable_names =          [ 'visc_type', 'plastic_type','softening_type' , 'density_type' ]
	variable_types =          [ 'int'      , 'int'         , 'int'           , 'int'          ]
	variable_extents        = [ 1          , 1             , 1               , 1              ]
	variable_textural_names = [ 'visc_type', 'plastic_type', 'softening_type', 'density_type' ]
	
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_DensityConst():
	ClassName      = 'MaterialConst_DensityConst'
	ClassNameShort = 'DensityConst'
	variable_names =          [ 'density']
	variable_types =          [ 'double' ]
	variable_extents        = [ 1        ]
	variable_textural_names = [ 'density']
        
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_DensityBoussinesq():
	ClassName      = 'MaterialConst_DensityBoussinesq'
	ClassNameShort = 'DensityBoussinesq'
	variable_names =          [ 'density','alpha'           ,    'beta'       ]
	variable_types =          [ 'double' ,'double'          ,    'double'     ]
	variable_extents        = [ 1        ,      1           ,       1         ]
	variable_textural_names = [ 'density','thermalexpension','compressibility']
	
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_ViscosityArrh():
	ClassName      = 'MaterialConst_ViscosityArrh'
	ClassNameShort = 'ViscosityArrh'
	variable_names =          [ 'preexpA','Ascale' ,'entalpy' , 'Vmol'   ,'nexp'    ,'Tref'    ,'Eta_scale','P_scale']       
	variable_types =          [ 'double' ,'double' , 'double' , 'double' , 'double' , 'double' ,'double'   ,'double' ]
	variable_extents        = [ 1        ,      1  ,       1  , 1        , 1        ,  1       , 1         ,    1    ]
	variable_textural_names =  [ 'preexpA','Ascale' ,'entalpy' , 'Vmol'  ,'nexp'    ,'Tref'    ,'Eta_scale','P_scale']       
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_ViscosityFK():
	ClassName      = 'MaterialConst_ViscosityFK'
	ClassNameShort = 'ViscosityFK'
	variable_names =          [ 'eta0'   , 'theta'  ]
	variable_types =          [ 'double' , 'double' ]
	variable_extents        = [ 1        ,    1     ]
	variable_textural_names = [ 'eta0'   ,'theta'   ]
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_SoftLin():
	ClassName      = 'MaterialConst_SoftLin'
	ClassNameShort = 'SoftLin'
	variable_names =          [ 'eps_min','eps_max'    ]
	variable_types =          [ 'double' ,'double'     ]
	variable_extents        = [ 1        ,      1      ]
	variable_textural_names = [ 'eps_min','eps_max'    ]
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )

def Generate_pTatin_MaterialConst_SoftExpo():
	ClassName      = 'MaterialConst_SoftExpo'
	ClassNameShort = 'SoftExpo'
	variable_names =          [ 'eps_min','eps_fold'   ]
	variable_types =          [ 'double' ,'double'     ]
	variable_extents        = [ 1        ,      1      ]
	variable_textural_names = [ 'eps_min','eps_fold'   ]
	PARTICLE_CLASS_GENERATOR( ClassName, ClassNameShort, variable_names, variable_types, variable_extents, variable_textural_names )


# Call all functions to generate all data types

## quadrature point fields ##
# Stokes
Generate_pTatin3d_QuadraturePointVolumeCoefficientStokes()
Generate_pTatin3d_QuadraturePointSurfaceCoefficientStokes()

# Energy
Generate_pTatin_QuadraturePointVolumeCoefficientEnergy()

## material point fields ##
Generate_pTatin3d_MaterialPointStandard()
Generate_pTatin_MaterialPointStokesData()
Generate_pTatin_MaterialPointStokesPlastic()
Generate_pTatin_MaterialPointEnergy()
Generate_pTatin_MaterialPointViscoElasticity()

## material constants ##
Generate_pTatin_MaterialConst_ViscosityConst()
Generate_pTatin_MaterialConst_ViscosityZ()
Generate_pTatin_MaterialConst_ViscosityArrh()
Generate_pTatin_MaterialConst_ViscosityFK()

Generate_pTatin_MaterialConst_DensityConst()
Generate_pTatin_MaterialConst_DensityBoussinesq()
Generate_pTatin_MaterialConst_PlasticMises()
Generate_pTatin_MaterialConst_PlasticDP()
Generate_pTatin_MaterialConst_SoftLin()
Generate_pTatin_MaterialConst_SoftExpo()
Generate_pTatin_MaterialConst_MaterialType()

# quad for surface processes
Generate_pTatin_QuadraturePointVolumeCoefficientSPM_a()
Generate_pTatin_QuadraturePointVolumeCoefficientSPM_standard_fields()
Generate_pTatin_QuadraturePointVolumeCoefficientSPM_nonlinear_diffusivity()

