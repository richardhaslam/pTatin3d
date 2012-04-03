## ===================================================
#
#    Affine transformations for 
#    2d (quads) and 3d (hex) elements
#
## ===================================================


from sympy import *
from sympy import var, Plot
from sympy import collect, sympify, Wild

# define variables we will differentiate with respect to
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

v1 = Symbol('v1')
v2 = Symbol('v2')
v3 = Symbol('v3')
v4 = Symbol('v4')

v5 = Symbol('v5')
v6 = Symbol('v6')
v7 = Symbol('v7')
v8 = Symbol('v8')


b1 = 0.25 * (1-x)*(1-y)
b2 = 0.25 * (1+x)*(1-y)
b3 = 0.25 * (1+x)*(1+y)
b4 = 0.25 * (1-x)*(1+y)

map = v1*b1 + v2*b2 + v3*b3 + v4*b4
em = expand(map)
#
#d = collect(em, [1,x,y], evaluate=False)
#print d[x]


print '2D affine map'

a00 = em.subs([(x,0),(y,0)]).evalf()
print 'a00 = ', a00

a10 = em.subs([(y,0)]).evalf()
a10 = a10 - a00
print 'a10 = ', a10.subs(x,1).evalf()

a01 = em.subs([(x,0)]).evalf()
a01 = a01 - a00
print 'a01 = ', a01.subs(y,1).evalf()




t1 = 0.125 * (1-x)*(1-y)*(1-z)
t2 = 0.125 * (1+x)*(1-y)*(1-z)
t3 = 0.125 * (1+x)*(1+y)*(1-z)
t4 = 0.125 * (1-x)*(1+y)*(1-z)

t5 = 0.125 * (1-x)*(1-y)*(1+z)
t6 = 0.125 * (1+x)*(1-y)*(1+z)
t7 = 0.125 * (1+x)*(1+y)*(1+z)
t8 = 0.125 * (1-x)*(1+y)*(1+z)

tmap = v1*t1 + v2*t2 + v3*t3 + v4*t4  +  v5*t5 + v6*t6 + v7*t7 + v8*t8
em = expand(tmap)


print '3D affine map'

a000 = em.subs([(x,0),(y,0),(z,0)]).evalf()
print 'a000 = ', a000


a100 = em.subs([(y,0),(z,0)]).evalf()
a100 = a100 - a000
print 'a100 = ', a100.subs(x,1).evalf()


a010 = em.subs([(x,0),(z,0)]).evalf()
a010 = a010 - a000
print 'a010 = ', a010.subs(y,1).evalf()


a001 = em.subs([(x,0),(y,0)]).evalf()
a001 = a001 - a000
print 'a001 = ', a001.subs(z,1).evalf()








