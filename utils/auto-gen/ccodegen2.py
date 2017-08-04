

class include:
  def __init__(self,filename,sys=False):
    self.filename=filename
    self.sys=True if sys==True else False
  
  def __str__(self):
    if self.sys == True:
      return '#include <%s>'%self.filename
    else:
      return '#include "%s"'%self.filename

class function:
  def __init__(self,name,argList,returnType):
    self.name=name
    self.argList = argList
    self.returnType = returnType
    self.stack = None
    self.buffer = ''

  def open(self):
    allargs = ','.join(self.argList)
    line = self.returnType + ' ' + self.name + '(' + str(allargs) + ')\n'

    self.stack = stackframe(None)
    line += self.stack.open()
    
    if self.returnType != 'void':
      line += self.private_statement([self.returnType + ' fReturnValue;'])
      line += '\n'
    
    self.buffer += line
    return(line)

  def petscopen(self):
    allargs = ','.join(self.argList)
    line = self.returnType + ' ' + self.name + '(' + str(allargs) + ')\n'
    
    self.stack = stackframe(None)
    line += self.stack.open()
    
    line += self.private_statement(['PetscErrorCode ierr;'])
    line += '\n'
    
    self.buffer += line
    return(line)

  # Method is identical to statement(), except this one does not buffer the result
  def private_statement(self,s):
    spacer = ' '*self.stack.indent
    line = ''
    for c in s:
      line += spacer + c + '\n'
    return(line)

  def statement(self,s):
    spacer = ' '*self.stack.indent
    #code = ';\n'.join(s)
    #code += code + ';\n'
    line = ''
    for c in s:
      line += spacer + c + '\n'
    self.buffer += line
    return(line)

  def close(self):
    line = ''
    if self.returnType == 'void':
      doNothing = True
      #line = self.private_statement(['return;'])
    else:
      line += self.private_statement(['return(fReturnValue);'])
    
    line += self.stack.close()
    self.buffer += line
    return(line)

  def petscclose(self):
    line = self.private_statement(['PetscFunctionReturn(0);'])
    line += self.stack.close()
    self.buffer += line
    return(line)

  def prototype(self):
    allargs = ','.join(self.argList)
    line = self.returnType + ' ' + self.name + '(' + str(allargs) + ');'
    return(line)

  def forBegin(self,index,range):
    sf = stackframe(self.stack)
    sf.prev = self.stack
    self.stack = sf
    
    cmd = 'for (' + index + '=' + str(range[0]) + '; ' + index + '<' + str(range[1]) + '; ' + index + '++)'
    code = self.private_statement([cmd])
    code += sf.open()
    self.buffer += code
    return(code)

  def forEnd(self,):
    currstack = self.stack
    code = currstack.close()
    self.stack = currstack.parent
    self.buffer += code
    return(code)

  def ifBegin(self,condition):
    sf = stackframe(self.stack)
    sf.prev = self.stack
    self.stack = sf
    
    cmd = 'if (' + str(condition) + ')'
    code = self.private_statement([cmd])
    code += sf.open()
    self.buffer += code
    return(code)

  def elseIfBegin(self,condition):
    currstack = self.stack
    code = currstack.close()
    self.stack = currstack.parent
    
    sf = stackframe(self.stack)
    sf.prev = self.stack
    self.stack = sf
    
    cmd = 'else if (' + str(condition) + ')'
    code += self.private_statement([cmd])
    code += sf.open()
    self.buffer += code
    return(code)

  def elseBegin(self):
    currstack = self.stack
    code = currstack.close()
    self.stack = currstack.parent
    
    sf = stackframe(self.stack)
    sf.prev = self.stack
    self.stack = sf
    
    cmd = 'else'
    code += self.private_statement([cmd])
    code += sf.open()
    self.buffer += code
    return(code)

  def ifEnd(self):
    currstack = self.stack
    code = currstack.close()
    self.stack = currstack.parent
    self.buffer += code
    return(code)

  def write(self):
    return(self.buffer)

class stackframe:
  def __init__(self,parent):
    self.parent = None
    self.indent = 0
    if parent != None:
      self.parent = parent
      self.indent = parent.indent
  
  def open(self):
    spacer = ' '*self.indent
    line = spacer + '{\n'
    self.indent += 2
    return(line)

  # note sure this method is required
  def statement(self,s):
    spacer = ' '*self.indent
    return(spacer + s + '\n')

  def close(self):
    self.indent -= 2
    spacer = ' '*self.indent
    line = spacer + '}\n'
    return(line)

class struct:
  def __init__(self,name,typeDef):
    self.indent = 0
    self.name = name
    self.typeDef = typeDef
  
  def open(self):
    spacer = ' '*self.indent
    if self.typeDef == False:
      line = spacer + 'struct ' + self.name + ' {\n'
    else:
      line = spacer + 'typedef struct {\n'
    self.indent += 2
    return(line)
  
  # note sure this method is required
  def statement(self,s):
    spacer = ' '*self.indent
    return(spacer + s + '\n')
  
  def close(self):
    self.indent -= 2
    spacer = ' '*self.indent
    if self.typeDef == False:
      line = spacer + '};\n'
    else:
      line = spacer + '} ' + self.name + ';\n'
    
    return(line)

