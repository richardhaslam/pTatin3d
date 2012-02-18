#
#
# Filter to define these operators
#
# [y0] = [ A11 A12 ] [u,v,w]
# [y1]   [ A21 0   ] [p]
#
# y0 = A11.[u,v,w]
# y0 = A12.[p]
# y1 = A21.[u,v,w]
#
# y0 = A11_xx.[u]
# y0 = A11_yy.[v]
#
# y0_x = A11.[0,v,w]
#
#

from sympy import *
from sympy import Function, dsolve, Eq, Derivative
from sympy.matrices import Matrix, zeros


def count_operations(var):
	ops = 0

	# count instances of +,-,*,/
	ops = ops + str(var).count('+')
	ops = ops + str(var).count('-')
	ops = ops + str(var).count('*')
	ops = ops + str(var).count('/')

	return ops


def delay_evaluation(var,name):
	return name

def ccode_array(var,nrow,ncol,outputVarName,add):

	if nrow != 1:
		if ncol != 1:
			for i in range(0,nrow):
				for j in range(0,ncol):
					if add == True:
						print outputVarName + '[' + str(i)+']['+str(j)+'] += ' + ccode(var[i,j]) + ';'
					else:
						print outputVarName + '[' + str(i)+']['+str(j)+'] = ' + ccode(var[i,j]) + ';'

	if nrow == 1:
		if ncol != 1:
			for i in range(0,ncol):
				if add == True:
					print outputVarName + '[' + str(i)+'] += ' + ccode(var[0,i]) + ';'
				else:
					print outputVarName + '[' + str(i)+'] = ' + ccode(var[0,i]) + ';'
				

	if ncol == 1:
		if nrow != 1:
			for i in range(0,nrow):
				if add == True:
					print outputVarName + '[' + str(i)+'] += ' + ccode(var[i,0]) + ';'
				else:
					print outputVarName + '[' + str(i)+'] = ' + ccode(var[i,0]) + ';'

	if nrow == 1:
		if ncol == 1:
			if add == True:
				print outputVarName + ' += ' + ccode(var[0,0]) + ';'
			else:
				print outputVarName + ' = ' + ccode(var[0,0]) + ';'

	return

# =========================================
def Stokes3dMixedFEM_MatMult_gp(nu,np,X,Y):

	eta_gp = Symbol('eta_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets

	for i in range(0,nsd*nu+np):
		Y[i,0] = 0

	# compute velocity and velocity gradients
	Uxdof = zeros((nu,1))
	Uydof = zeros((nu,1))
	Uzdof = zeros((nu,1))
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	UUdof = zeros((nsd*nu,1))
	for i in range(0,nu):
		UUdof[nsd*i  ,0] = X[i  ]
		UUdof[nsd*i+1,0] = X[i+nu]
		UUdof[nsd*i+2,0] = X[i+nu+nu]

		Uxdof[i] = UUdof[nsd*i  ,0]
		Uydof[i] = UUdof[nsd*i+1,0]
		Uzdof[i] = UUdof[nsd*i+2,0]

		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

	eval_ux_gp = Nu * Uxdof
	eval_uy_gp = Nu * Uydof
	eval_uz_gp = Nu * Uzdof

	print '// velocity(x) at gauss point //'
	ccode_array(eval_ux_gp,1,1,'u_gp',False)
	print '// velocity(y) at gauss point //'
	ccode_array(eval_uy_gp,1,1,'v_gp',False)
	print '// velocity(z) at gauss point //'
	ccode_array(eval_uz_gp,1,1,'w_gp',False)


	dudx_gp = dNudx.transpose() * Uxdof
	dudy_gp = dNudy.transpose() * Uxdof
	dudz_gp = dNudz.transpose() * Uxdof

	dvdx_gp = dNudx.transpose() * Uydof
	dvdy_gp = dNudy.transpose() * Uydof
	dvdz_gp = dNudz.transpose() * Uydof

	dwdx_gp = dNudx.transpose() * Uzdof
	dwdy_gp = dNudy.transpose() * Uzdof
	dwdz_gp = dNudz.transpose() * Uzdof

		
	# compute pressure		
	Pdof = zeros((np,1))
	Np   = zeros((1,np))

	for i in range(0,np):
		Pdof[i,0] = X[nsd*nu+i]
		Np[0,i] = Symbol('Np['+str(i)+']')

	eval_p_gp  = Np * Pdof
	print '// pressure at gauss point //'
	ccode_array(eval_p_gp,1,1,'p_gp',False)
	
	# compute strain-rate
	eval_E_gp = zeros((ntens,1))
	eval_E_gp[0,0] = 0.5 * ( dudx_gp[0,0] + dudx_gp[0,0] )
	eval_E_gp[1,0] = 0.5 * ( dvdy_gp[0,0] + dvdy_gp[0,0] )
	eval_E_gp[2,0] = 0.5 * ( dwdz_gp[0,0] + dwdz_gp[0,0] )
	
	eval_E_gp[3,0] = ( dudy_gp[0,0] + dvdx_gp[0,0] )
	eval_E_gp[4,0] = ( dudz_gp[0,0] + dwdx_gp[0,0] )
	eval_E_gp[5,0] = ( dvdz_gp[0,0] + dwdy_gp[0,0] )


	print '// strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'
	ccode_array(eval_E_gp,ntens,1,'strain_rate_gp',False)

	E_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_E_gp[d,0] != 0:
			E_gp[d,0] = Symbol('strain_rate_gp['+str(d)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	# compute deviatoric stress
	eval_tau_gp = D * E_gp 
	print '// deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //'
	ccode_array(eval_tau_gp,ntens,1,'tau_gp',False)

	tau_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_tau_gp[d,0] != 0:
			tau_gp[d,0] = Symbol('tau_gp['+str(d)+']')

	## Bt.D.B.u
	B = zeros((ntens,nsd*nu))
	for i in range(0,nu):
		B[0,nsd*i+0] = Symbol('dNudx['+str(i)+']')
		B[1,nsd*i+1] = Symbol('dNudy['+str(i)+']')
		B[2,nsd*i+2] = Symbol('dNudz['+str(i)+']')

		# exy = du/dy + dv/dx
		B[3,nsd*i+0] = Symbol('dNudy['+str(i)+']')
		B[3,nsd*i+1] = Symbol('dNudx['+str(i)+']')

		# exz = du/dz + dw/dx
		B[4,nsd*i+0] = Symbol('dNudz['+str(i)+']')
		B[4,nsd*i+2] = Symbol('dNudx['+str(i)+']')

		# eyz = dv/dz + dw/dy
		B[5,nsd*i+1] = Symbol('dNudz['+str(i)+']')
		B[5,nsd*i+2] = Symbol('dNudy['+str(i)+']')


	Bt_i_tau_gp = B.transpose() * tau_gp
	print '// y = A11.u at gauss point //'
	ccode_array(Bt_i_tau_gp,nsd*nu,1,'y0_uu',True)


	## -Bt.p.m
	m = zeros((ntens,1))
	m[0,0] = 1
	m[1,0] = 1
	m[2,0] = 1
	
	m[3,0] = 0
	m[4,0] = 0
	m[5,0] = 0

	p_gp = 0
	if eval_p_gp[0,0] != 0:
		p_gp = Symbol('p_gp')

	p_m_gp = p_gp * m

	grad_pressure_gp = -B.transpose() * p_m_gp
	print '// y = A12.p at gauss point //'
	ccode_array(grad_pressure_gp,nsd*nu,1,'y0_up',True)


	## Need to transpose this operator; -Bt.m.Np.p^h
	## -Np^T.m^T.B.u
	## -Np_i^T.m^T.(B.u)_j
	eval_div_gp = m.transpose() * E_gp
	print '// divergence at gauss point //'
	ccode_array( eval_div_gp,1,1, 'div_gp',True)

	div_gp = 0
	if eval_div_gp[0,0] != 0:
		div_gp = Symbol('div_gp')

	Np_div_u = -Np.transpose() * div_gp
	
	print '// y = A21.u at gauss point //'
	ccode_array(Np_div_u,np,1,'y1_pu',True)

	# combine results into one vector
	for i in range(nu):
		Y[nsd*i  ] = Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0]
		Y[nsd*i+1] = Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0]
		Y[nsd*i+2] = Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0]

	for i in range(np):
		Y[nsd*nu+i] = Np_div_u[i,0]

# =========================================
def Stokes3dStabFEM_MatMult_gp(nu,np,X,Y):

	eta_gp = Symbol('eta_gp')
	eta_avg_gp = Symbol('eta_avg_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets

	for i in range(0,nsd*nu+np):
		Y[i,0] = 0

	# compute velocity and velocity gradients
	Uxdof = zeros((nu,1))
	Uydof = zeros((nu,1))
	Uzdof = zeros((nu,1))
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	UUdof = zeros((nsd*nu,1))
	for i in range(0,nu):
		UUdof[nsd*i  ,0] = X[i  ]
		UUdof[nsd*i+1,0] = X[i+nu]
		UUdof[nsd*i+2,0] = X[i+nu+nu]

		Uxdof[i] = UUdof[nsd*i  ,0]
		Uydof[i] = UUdof[nsd*i+1,0]
		Uzdof[i] = UUdof[nsd*i+2,0]

		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

	eval_ux_gp = Nu * Uxdof
	eval_uy_gp = Nu * Uydof
	eval_uz_gp = Nu * Uzdof

	print '// velocity(x) at gauss point //'
	ccode_array(eval_ux_gp,1,1,'u_gp',False)
	print '// velocity(y) at gauss point //'
	ccode_array(eval_uy_gp,1,1,'v_gp',False)
	print '// velocity(z) at gauss point //'
	ccode_array(eval_uz_gp,1,1,'w_gp',False)


	dudx_gp = dNudx.transpose() * Uxdof
	dudy_gp = dNudy.transpose() * Uxdof
	dudz_gp = dNudz.transpose() * Uxdof

	dvdx_gp = dNudx.transpose() * Uydof
	dvdy_gp = dNudy.transpose() * Uydof
	dvdz_gp = dNudz.transpose() * Uydof

	dwdx_gp = dNudx.transpose() * Uzdof
	dwdy_gp = dNudy.transpose() * Uzdof
	dwdz_gp = dNudz.transpose() * Uzdof

		
	# compute pressure		
	Pdof = zeros((np,1))
	Np   = zeros((1,np))

	for i in range(0,np):
		Pdof[i,0] = X[nsd*nu+i]
		Np[0,i] = Symbol('Np['+str(i)+']')

	eval_p_gp  = Np * Pdof
	print '// pressure at gauss point //'
	ccode_array(eval_p_gp,1,1,'p_gp',False)
	
	# compute strain-rate
	eval_E_gp = zeros((ntens,1))
	eval_E_gp[0,0] = 0.5 * ( dudx_gp[0,0] + dudx_gp[0,0] )
	eval_E_gp[1,0] = 0.5 * ( dvdy_gp[0,0] + dvdy_gp[0,0] )
	eval_E_gp[2,0] = 0.5 * ( dwdz_gp[0,0] + dwdz_gp[0,0] )
	
	eval_E_gp[3,0] = ( dudy_gp[0,0] + dvdx_gp[0,0] )
	eval_E_gp[4,0] = ( dudz_gp[0,0] + dwdx_gp[0,0] )
	eval_E_gp[5,0] = ( dvdz_gp[0,0] + dwdy_gp[0,0] )


	print '// strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'
	ccode_array(eval_E_gp,ntens,1,'strain_rate_gp',False)

	E_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_E_gp[d,0] != 0:
			E_gp[d,0] = Symbol('strain_rate_gp['+str(d)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	# compute deviatoric stress
	eval_tau_gp = D * E_gp 
	print '// deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //'
	ccode_array(eval_tau_gp,ntens,1,'tau_gp',False)

	tau_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_tau_gp[d,0] != 0:
			tau_gp[d,0] = Symbol('tau_gp['+str(d)+']')

	## Bt.D.B.u
	B = zeros((ntens,nsd*nu))
	for i in range(0,nu):
		B[0,nsd*i+0] = Symbol('dNudx['+str(i)+']')
		B[1,nsd*i+1] = Symbol('dNudy['+str(i)+']')
		B[2,nsd*i+2] = Symbol('dNudz['+str(i)+']')

		# exy = du/dy + dv/dx
		B[3,nsd*i+0] = Symbol('dNudy['+str(i)+']')
		B[3,nsd*i+1] = Symbol('dNudx['+str(i)+']')

		# exz = du/dz + dw/dx
		B[4,nsd*i+0] = Symbol('dNudz['+str(i)+']')
		B[4,nsd*i+2] = Symbol('dNudx['+str(i)+']')

		# eyz = dv/dz + dw/dy
		B[5,nsd*i+1] = Symbol('dNudz['+str(i)+']')
		B[5,nsd*i+2] = Symbol('dNudy['+str(i)+']')


	Bt_i_tau_gp = B.transpose() * tau_gp
	print '// y = A11.u at gauss point //'
	ccode_array(Bt_i_tau_gp,nsd*nu,1,'y0_uu',True)


	## -Bt.p.m
	m = zeros((ntens,1))
	m[0,0] = 1
	m[1,0] = 1
	m[2,0] = 1
	
	m[3,0] = 0
	m[4,0] = 0
	m[5,0] = 0

	p_gp = 0
	if eval_p_gp[0,0] != 0:
		p_gp = Symbol('p_gp')

	p_m_gp = p_gp * m

	grad_pressure_gp = -B.transpose() * p_m_gp
	print '// y = A12.p at gauss point //'
	ccode_array(grad_pressure_gp,nsd*nu,1,'y0_up',True)


	## Need to transpose this operator; -Bt.m.Np.p^h
	## -Np^T.m^T.B.u
	## -Np_i^T.m^T.(B.u)_j
	eval_div_gp = m.transpose() * E_gp
	print '// divergence at gauss point //'
	ccode_array( eval_div_gp,1,1, 'div_gp',True)

	div_gp = 0
	if eval_div_gp[0,0] != 0:
		div_gp = Symbol('div_gp')

	Np_div_u = -Np.transpose() * div_gp
	
	print '// y = A21.u at gauss point //'
	ccode_array(Np_div_u,np,1,'y1_pu',True)


	print '// shifted pressure at gauss point //'
	
	ones = zeros((np,1))
	for d in range(0,np):
		ones[d,0] = 1.0

	eval_sum_p = Pdof.transpose() * ones
	ccode_array( eval_sum_p,1,1,'sum_p',False)

	sum_p = 0
	if eval_sum_p[0,0] != 0:
		sum_p = Symbol('sum_p')

	eval_shifted_p_gp = (1.0/eta_avg_gp)*( p_gp - (1.0/18.0) * sum_p )
	#ccode_array( shifted_p_gp,1,1,'shifted_p_gp',False)
	print 'shifted_p_gp = ' + str(eval_shifted_p_gp) + ';'

	shifted_p_gp = 0
	if eval_shifted_p_gp != 0:
		shifted_p_gp = Symbol('shifted_p_gp')

	stab = zeros((np,1))
	stab = - Np.transpose() * shifted_p_gp
	ccode_array(stab,np,1,'y1_pp',True)



	# combine results into one vector
	for i in range(nu):
		Y[nsd*i  ] = Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0]
		Y[nsd*i+1] = Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0]
		Y[nsd*i+2] = Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0]

	for i in range(np):
		Y[nsd*nu+i] = Np_div_u[i,0] + stab[i,0]


def Stokes3dMixedFEM_MatMult_gp_with_filter(funcname,nu,np,X,vfields,pfield):

	eta_gp = Symbol('eta_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets


	print 'inline void ' + funcname + '(const double FAC,const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])'
	print '{'
	print '  const int nsd = '+str(nsd)+';'
	print '  const int ntens = '+str(ntens)+';'
	print '  double    p_gp;'
	print '  double    strain_rate_gp['+str(ntens)+'];'
	print '  double    tau_gp['+str(ntens)+'];'
	print '  double    div_gp;'
#	print '  double Nu['+str(nu)+'];'
#	print '  double Np['+str(np)+'];'
#	print '  double dNudx['+str(nu)+'];'
#	print '  double dNudy['+str(nu)+'];'
	print '\n'

	Y = zeros((nsd*nu+np,1))
	total_ops = 0

	# compute velocity and velocity gradients
	Uxdof = zeros((nu,1))
	Uydof = zeros((nu,1))
	Uzdof = zeros((nu,1))
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	UUdof = zeros((nsd*nu,1))
	for i in range(0,nu):
		UUdof[nsd*i  ,0] = X[i  ]
		UUdof[nsd*i+1,0] = X[i+nu]
		UUdof[nsd*i+2,0] = X[i+nu+nu]

		Uxdof[i] = UUdof[nsd*i  ,0]
		Uydof[i] = UUdof[nsd*i+1,0]
		Uzdof[i] = UUdof[nsd*i+2,0]

		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

	eval_ux_gp = Nu * Uxdof
	eval_uy_gp = Nu * Uydof
	eval_uz_gp = Nu * Uzdof

	#print '  // velocity(x) at gauss point //'
	#ccode_array(eval_ux_gp,1,1,'  u_gp')
	#print '  // velocity(y) at gauss point //'
	#ccode_array(eval_uy_gp,1,1,'  v_gp')
	#print '  // velocity(z) at gauss point //'
	#ccode_array(eval_uz_gp,1,1,'  w_gp')

	#ops = count_operations(eval_ux_gp)
	#ops = ops + count_operations(eval_uy_gp)
	#ops = ops + count_operations(eval_uz_gp)
	#total_ops = total_ops + ops

	dudx_gp = dNudx.transpose() * Uxdof
	dudy_gp = dNudy.transpose() * Uxdof
	dudz_gp = dNudz.transpose() * Uxdof

	dvdx_gp = dNudx.transpose() * Uydof
	dvdy_gp = dNudy.transpose() * Uydof
	dvdz_gp = dNudz.transpose() * Uydof

	dwdx_gp = dNudx.transpose() * Uzdof
	dwdy_gp = dNudy.transpose() * Uzdof
	dwdz_gp = dNudz.transpose() * Uzdof

		
	# compute pressure		
	Pdof = zeros((np,1))
	Np   = zeros((1,np))

	for i in range(0,np):
		Pdof[i,0] = X[nsd*nu+i]
		Np[0,i] = Symbol('Np['+str(i)+']')

	eval_p_gp  = Np * Pdof
	print '  // pressure at gauss point //'
	ccode_array(eval_p_gp,1,1,'  p_gp',False)

	ops = count_operations(eval_p_gp)
	total_ops = total_ops + ops
	
	# compute strain-rate
	eval_E_gp = zeros((ntens,1))
	eval_E_gp[0,0] = 0.5 * ( dudx_gp[0,0] + dudx_gp[0,0] )
	eval_E_gp[1,0] = 0.5 * ( dvdy_gp[0,0] + dvdy_gp[0,0] )
	eval_E_gp[2,0] = 0.5 * ( dwdz_gp[0,0] + dwdz_gp[0,0] )
	
	eval_E_gp[3,0] = ( dudy_gp[0,0] + dvdx_gp[0,0] )
	eval_E_gp[4,0] = ( dudz_gp[0,0] + dwdx_gp[0,0] )
	eval_E_gp[5,0] = ( dvdz_gp[0,0] + dwdy_gp[0,0] )


	print '  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'
	ccode_array(eval_E_gp,ntens,1,'  strain_rate_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_E_gp[d,0])

	total_ops = total_ops + ops

	E_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_E_gp[d,0] != 0:
			E_gp[d,0] = Symbol('strain_rate_gp['+str(d)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	# compute deviatoric stress
	eval_tau_gp = D * E_gp 
	print '  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //'
	ccode_array(eval_tau_gp,ntens,1,'  tau_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_tau_gp[d,0])

	total_ops = total_ops + ops

	tau_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_tau_gp[d,0] != 0:
			tau_gp[d,0] = Symbol('tau_gp['+str(d)+']')

	## Bt.D.B.u
	B = zeros((ntens,nsd*nu))
	for i in range(0,nu):
		B[0,nsd*i+0] = Symbol('dNudx['+str(i)+']')
		B[1,nsd*i+1] = Symbol('dNudy['+str(i)+']')
		B[2,nsd*i+2] = Symbol('dNudz['+str(i)+']')

		# exy = du/dy + dv/dx
		B[3,nsd*i+0] = Symbol('dNudy['+str(i)+']')
		B[3,nsd*i+1] = Symbol('dNudx['+str(i)+']')

		# exz = du/dz + dw/dx
		B[4,nsd*i+0] = Symbol('dNudz['+str(i)+']')
		B[4,nsd*i+2] = Symbol('dNudx['+str(i)+']')

		# eyz = dv/dz + dw/dy
		B[5,nsd*i+1] = Symbol('dNudz['+str(i)+']')
		B[5,nsd*i+2] = Symbol('dNudy['+str(i)+']')


	Bt_i_tau_gp = B.transpose() * tau_gp
	print '  // y = A11.u at gauss point //'
	#ccode_array(Bt_i_tau_gp,nsd*nu,1,'y0_uu')

	## -Bt.p.m
	m = zeros((ntens,1))
	m[0,0] = 1
	m[1,0] = 1
	m[2,0] = 1
	
	m[3,0] = 0
	m[4,0] = 0
	m[5,0] = 0

	p_gp = 0
	if eval_p_gp[0,0] != 0:
		p_gp = Symbol('p_gp')

	p_m_gp = p_gp * m

	grad_pressure_gp = -B.transpose() * p_m_gp
	print '  // y = A12.p at gauss point //'
	#ccode_array(grad_pressure_gp,nsd*nu,1,'y0_up')


	## Need to transpose this operator; -Bt.m.Np.p^h
	## -Np^T.m^T.B.u
	## -Np_i^T.m^T.(B.u)_j
	eval_div_gp = m.transpose() * E_gp
	print '  // divergence at gauss point //'
	ccode_array( eval_div_gp,1,1, '  div_gp', False )
	
	ops = count_operations(eval_div_gp)
	total_ops = total_ops + ops

	div_gp = 0
	if eval_div_gp[0,0] != 0:
		div_gp = Symbol('div_gp')

	Np_div_u = -Np.transpose() * div_gp
	
	print '  // y = A21.u at gauss point //'
	#ccode_array(Np_div_u,np,1,'y1_pu')

#	# combine results into one vector
#	for i in range(nu):
#		Y[nsd*i  ] = Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0]
#		Y[nsd*i+1] = Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0]
#		Y[nsd*i+2] = Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0]

#	for i in range(np):
#		Y[nsd*nu+i] = Np_div_u[i,0]

	# combine results into one vector
	FAC = Symbol('FAC')
	for i in range(nu):
		Y[nsd*i  ] = FAC*(Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0])
		Y[nsd*i+1] = FAC*(Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0])
		Y[nsd*i+2] = FAC*(Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0])

	for i in range(np):
		Y[nsd*nu+i] = FAC*Np_div_u[i,0]

	Yextract = extract_fields(Y,nsd,nu,vfields,np,pfield)
	lu = len(vfields)
	lp = len(pfield)
	
	ops = 0
	for d in range(0,lu*nu+lp*np):
		ops = ops + count_operations(Yextract[d,0])
	total_ops = total_ops + ops

	# add += 
	total_ops = total_ops + (lu*nu+lp*np)


	ccode_array(Yextract,lu*nu+lp*np,1,'  Y',True)
	print '  \n// total operations = ' + str(total_ops)

	print '}'


def Stokes3dMixedFEM_diag_MatMult_gp_with_filter(funcname,nu,np,vfields,pfield):

	eta_gp = Symbol('eta_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets


	print 'inline void ' + funcname + '(const double FAC,const double eta_gp,const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])'
	print '{'
	print '  const int nsd = '+str(nsd)+';'
	print '  const int ntens = '+str(ntens)+';'
	print '\n'

	Y = zeros((nsd*nu+np,1))
	total_ops = 0

	# compute velocity and velocity gradients
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	for i in range(0,nu):
		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

		
	# compute pressure		
	Np   = zeros((1,np))

	for i in range(0,np):
		Np[0,i] = Symbol('Np['+str(i)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	
	# compute strain-rate operator
	eval_B1_gp = zeros((ntens,nu))
	eval_B2_gp = zeros((ntens,nu))
	eval_B3_gp = zeros((ntens,nu))

	for i in range(0,nu):
		eval_B1_gp[0,i] = dNudx[i,0]
		eval_B1_gp[3,i] = dNudy[i,0]
		eval_B1_gp[4,i] = dNudz[i,0]

		eval_B2_gp[1,i] = dNudy[i,0]
		eval_B2_gp[3,i] = dNudx[i,0]
		eval_B2_gp[5,i] = dNudz[i,0]

		eval_B3_gp[2,i] = dNudz[i,0]
		eval_B3_gp[4,i] = dNudx[i,0]
		eval_B3_gp[5,i] = dNudy[i,0]

#	ccode_array(eval_B1_gp,ntens,nu,'  B1_gp',False)
#	ccode_array(eval_B2_gp,ntens,nu,'  B1_gp',False)
#	ccode_array(eval_B3_gp,ntens,nu,'  B1_gp',False)


#	print '  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'

	Yu = zeros((nu,1))
	Yv = zeros((nu,1))
	Yw = zeros((nu,1))
	for ii in range(0,nu):
		a = 0
		b = 0
		c = 0
		for ss in range(0,ntens):
			for tt in range(0,ntens):
				a = a + eval_B1_gp[ss,ii] * D[ss,tt] * eval_B1_gp[tt,ii]
				b = b + eval_B2_gp[ss,ii] * D[ss,tt] * eval_B2_gp[tt,ii]
				c = c + eval_B3_gp[ss,ii] * D[ss,tt] * eval_B3_gp[tt,ii]
			
		
		Yu[ii] = a
		Yv[ii] = b
		Yw[ii] = c


#	Yu = eval_B1_gp.transpose() * D * eval_B1_gp
#	Yv = eval_B2_gp.transpose() * D * eval_B2_gp
#	Yw = eval_B3_gp.transpose() * D * eval_B3_gp

#	ccode_array(Yu,nu,1,'  Yu_gp',False)
#	ccode_array(Yv,nu,1,'  Yv_gp',False)
#	ccode_array(Yw,nu,1,'  Yw_gp',False)

	# combine results into one vector
	FAC = Symbol('FAC')
	for i in range(nu):
		Y[nsd*i  ] = FAC * Yu[i]
		Y[nsd*i+1] = FAC * Yv[i]
		Y[nsd*i+2] = FAC * Yw[i]

	for i in range(np):
		Y[nsd*nu+i] = FAC * 0

#	ccode_array(Y,3*nu+np,1,'  Y_gp',False)

	# prune components
	Yextract = extract_fields(Y,nsd,nu,vfields,np,pfield)
	lu = len(vfields)
	lp = len(pfield)
	
	ops = 0
	for d in range(0,lu*nu+lp*np):
		ops = ops + count_operations(Yextract[d,0])
	total_ops = total_ops + ops

	# add += 
	total_ops = total_ops + (lu*nu+lp*np)


	ccode_array(Yextract,lu*nu+lp*np,1,'  Y',True)
	print '  \n// total operations = ' + str(total_ops)




	print '}'

def Stokes3dStabFEM_MatMult_gp_with_filter(funcname,nu,np,X,vfields,pfield):

	eta_gp = Symbol('eta_gp')
	eta_avg_gp = Symbol('eta_avg_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets


	print 'inline void ' + funcname + '(const double eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])'
	print '{'
	print '  const int ntens = 6;'
	print '  double u_gp,v_gp,w_gp,p_gp;'
	print '  double strain_rate_gp[ntens];'
	print '  double tau_gp[ntens];'
	print '  double div_gp;'
	print '  int operations = 0;'
#	print '  double Nu['+str(nu)+'];'
#	print '  double Np['+str(np)+'];'
#	print '  double dNudx['+str(nu)+'];'
#	print '  double dNudy['+str(nu)+'];'
	print '\n'

	Y = zeros((nsd*nu+np,1))
	total_ops = 0

	# compute velocity and velocity gradients
	Uxdof = zeros((nu,1))
	Uydof = zeros((nu,1))
	Uzdof = zeros((nu,1))
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	UUdof = zeros((nsd*nu,1))
	for i in range(0,nu):
		UUdof[nsd*i  ,0] = X[i  ]
		UUdof[nsd*i+1,0] = X[i+nu]
		UUdof[nsd*i+2,0] = X[i+nu+nu]

		Uxdof[i] = UUdof[nsd*i  ,0]
		Uydof[i] = UUdof[nsd*i+1,0]
		Uzdof[i] = UUdof[nsd*i+2,0]

		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

	eval_ux_gp = Nu * Uxdof
	eval_uy_gp = Nu * Uydof
	eval_uz_gp = Nu * Uzdof

	#print '  // velocity(x) at gauss point //'
	#ccode_array(eval_ux_gp,1,1,'  u_gp')
	#print '  // velocity(y) at gauss point //'
	#ccode_array(eval_uy_gp,1,1,'  v_gp')
	#print '  // velocity(z) at gauss point //'
	#ccode_array(eval_uz_gp,1,1,'  w_gp')

	#ops = count_operations(eval_ux_gp)
	#ops = ops + count_operations(eval_uy_gp)
	#ops = ops + count_operations(eval_uz_gp)
	#total_ops = total_ops + ops

	dudx_gp = dNudx.transpose() * Uxdof
	dudy_gp = dNudy.transpose() * Uxdof
	dudz_gp = dNudz.transpose() * Uxdof

	dvdx_gp = dNudx.transpose() * Uydof
	dvdy_gp = dNudy.transpose() * Uydof
	dvdz_gp = dNudz.transpose() * Uydof

	dwdx_gp = dNudx.transpose() * Uzdof
	dwdy_gp = dNudy.transpose() * Uzdof
	dwdz_gp = dNudz.transpose() * Uzdof

		
	# compute pressure		
	Pdof = zeros((np,1))
	Np   = zeros((1,np))

	for i in range(0,np):
		Pdof[i,0] = X[nsd*nu+i]
		Np[0,i] = Symbol('Np['+str(i)+']')

	eval_p_gp  = Np * Pdof
	print '  // pressure at gauss point //'
	ccode_array(eval_p_gp,1,1,'  p_gp',False)

	ops = count_operations(eval_p_gp)
	total_ops = total_ops + ops
	
	# compute strain-rate
	eval_E_gp = zeros((ntens,1))
	eval_E_gp[0,0] = 0.5 * ( dudx_gp[0,0] + dudx_gp[0,0] )
	eval_E_gp[1,0] = 0.5 * ( dvdy_gp[0,0] + dvdy_gp[0,0] )
	eval_E_gp[2,0] = 0.5 * ( dwdz_gp[0,0] + dwdz_gp[0,0] )
	
	eval_E_gp[3,0] = ( dudy_gp[0,0] + dvdx_gp[0,0] )
	eval_E_gp[4,0] = ( dudz_gp[0,0] + dwdx_gp[0,0] )
	eval_E_gp[5,0] = ( dvdz_gp[0,0] + dwdy_gp[0,0] )


	print '  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'
	ccode_array(eval_E_gp,ntens,1,'  strain_rate_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_E_gp[d,0])

	total_ops = total_ops + ops

	E_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_E_gp[d,0] != 0:
			E_gp[d,0] = Symbol('strain_rate_gp['+str(d)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	# compute deviatoric stress
	eval_tau_gp = D * E_gp 
	print '  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //'
	ccode_array(eval_tau_gp,ntens,1,'  tau_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_tau_gp[d,0])

	total_ops = total_ops + ops

	tau_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_tau_gp[d,0] != 0:
			tau_gp[d,0] = Symbol('tau_gp['+str(d)+']')

	## Bt.D.B.u
	B = zeros((ntens,nsd*nu))
	for i in range(0,nu):
		B[0,nsd*i+0] = Symbol('dNudx['+str(i)+']')
		B[1,nsd*i+1] = Symbol('dNudy['+str(i)+']')
		B[2,nsd*i+2] = Symbol('dNudz['+str(i)+']')

		# exy = du/dy + dv/dx
		B[3,nsd*i+0] = Symbol('dNudy['+str(i)+']')
		B[3,nsd*i+1] = Symbol('dNudx['+str(i)+']')

		# exz = du/dz + dw/dx
		B[4,nsd*i+0] = Symbol('dNudz['+str(i)+']')
		B[4,nsd*i+2] = Symbol('dNudx['+str(i)+']')

		# eyz = dv/dz + dw/dy
		B[5,nsd*i+1] = Symbol('dNudz['+str(i)+']')
		B[5,nsd*i+2] = Symbol('dNudy['+str(i)+']')


	Bt_i_tau_gp = B.transpose() * tau_gp
	print '  // y = A11.u at gauss point //'
	#ccode_array(Bt_i_tau_gp,nsd*nu,1,'y0_uu')

	## -Bt.p.m
	m = zeros((ntens,1))
	m[0,0] = 1
	m[1,0] = 1
	m[2,0] = 1
	
	m[3,0] = 0
	m[4,0] = 0
	m[5,0] = 0

	p_gp = 0
	if eval_p_gp[0,0] != 0:
		p_gp = Symbol('p_gp')

	p_m_gp = p_gp * m

	grad_pressure_gp = -B.transpose() * p_m_gp
	print '  // y = A12.p at gauss point //'
	#ccode_array(grad_pressure_gp,nsd*nu,1,'y0_up')


	## Need to transpose this operator; -Bt.m.Np.p^h
	## -Np^T.m^T.B.u
	## -Np_i^T.m^T.(B.u)_j
	eval_div_gp = m.transpose() * E_gp
	print '  // divergence at gauss point //'
	ccode_array( eval_div_gp,1,1, '  div_gp', False )
	
	ops = count_operations(eval_div_gp)
	total_ops = total_ops + ops

	div_gp = 0
	if eval_div_gp[0,0] != 0:
		div_gp = Symbol('div_gp')

	Np_div_u = -Np.transpose() * div_gp
	
	print '  // y = A21.u at gauss point //'
	#ccode_array(Np_div_u,np,1,'y1_pu')


	print '  // shifted pressure at gauss point //'
	ones = zeros((np,1))
	for d in range(0,np):
		ones[d,0] = 1.0

	eval_sum_p = Pdof.transpose() * ones
	ccode_array( eval_sum_p,1,1,'  sum_p',False)
	ops = count_operations(eval_sum_p)
	total_ops = total_ops + ops

	sum_p = 0
	if eval_sum_p[0,0] != 0:
		sum_p = Symbol('sum_p')

	eval_shifted_p_gp = (1.0/eta_avg_gp)*( p_gp - (1.0/64.0) * sum_p )
	ops = count_operations(eval_shifted_p_gp)
	total_ops = total_ops + ops

	#ccode_array( shifted_p_gp,1,1,'shifted_p_gp',False)
	print '  shifted_p_gp = ' + str(eval_shifted_p_gp) + ';'

	shifted_p_gp = 0
	if eval_shifted_p_gp != 0:
		shifted_p_gp = Symbol('shifted_p_gp')

	stab = zeros((np,1))
	stab = - Np.transpose() * shifted_p_gp
#	ccode_array(stab,np,1,'  y1_pp',True)


	# combine results into one vector
	for i in range(nu):
		Y[nsd*i  ] = Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0]
		Y[nsd*i+1] = Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0]
		Y[nsd*i+2] = Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0]

	for i in range(np):
		Y[nsd*nu+i] = Np_div_u[i,0] + stab[i,0]

	Yextract = extract_fields(Y,nsd,nu,vfields,np,pfield)
	lu = len(vfields)
	lp = len(pfield)
	
	ops = 0
	for d in range(0,lu*nu+lp*np):
		ops = ops + count_operations(Yextract[d,0])
	total_ops = total_ops + ops

	# add += 
	total_ops = total_ops + (lu*nu+lp*np)


	ccode_array(Yextract,lu*nu+lp*np,1,'  Y',True)
	print '  \n// total operations = ' + str(total_ops)

	print '}'


##########################################################################################		
def Stokes3dStabFEM_MatMult_gp_with_filter_DEBUG(funcname,nu,np,X,vfields,pfield):

	eta_gp = Symbol('eta_gp')
	eta_avg_gp = Symbol('eta_avg_gp')
	nsd   = 3 # spatial dimensions 
	ntens = 6 # tensor componenets


	print 'inline void ' + funcname + '(const double FAC,const double eta_gp,const double avg_eta_gp,const double Ux[],const double Uy[],const double Uz[],const double P[],const double Nu[],const double dNudx[],const double dNudy[],const double dNudz[],const double Np[],double Y[])'
	print '{'
	
	print '  const int nsd = '+str(nsd)+';'
	print '  const int ntens = '+str(ntens)+';'
	print '  double    p_gp;'
	print '  double    strain_rate_gp['+str(ntens)+'];'
	print '  double    tau_gp['+str(ntens)+'];'
	print '  double    div_gp;'
	print '  int       operations = 0;'
	print '  double    div_gp;'
	print '  double    sum_p;'
	print '  double    scaled_p;'
	print '  double    one_on_eta_avg_gp;'
	print '  double    p_gp_on_eta_avg;'
	print '  double    scaled_p_on_eta_avg;'
#	print '  double Nu['+str(nu)+'];'
#	print '  double Np['+str(np)+'];'
#	print '  double dNudx['+str(nu)+'];'
#	print '  double dNudy['+str(nu)+'];'
	print '\n'

	Y = zeros((nsd*nu+np,1))
	total_ops = 0

	# compute velocity and velocity gradients
	Uxdof = zeros((nu,1))
	Uydof = zeros((nu,1))
	Uzdof = zeros((nu,1))
	Nu = zeros((1,nu))
	dNudx = zeros((nu,1))
	dNudy = zeros((nu,1))
	dNudz = zeros((nu,1))

	UUdof = zeros((nsd*nu,1))
	for i in range(0,nu):
		UUdof[nsd*i  ,0] = X[i  ]
		UUdof[nsd*i+1,0] = X[i+nu]
		UUdof[nsd*i+2,0] = X[i+nu+nu]

		Uxdof[i] = UUdof[nsd*i  ,0]
		Uydof[i] = UUdof[nsd*i+1,0]
		Uzdof[i] = UUdof[nsd*i+2,0]

		Nu[0,i] = Symbol('Nu['+str(i)+']')
		dNudx[i,0] = Symbol('dNudx['+str(i)+']')
		dNudy[i,0] = Symbol('dNudy['+str(i)+']')
		dNudz[i,0] = Symbol('dNudz['+str(i)+']')

	eval_ux_gp = Nu * Uxdof
	eval_uy_gp = Nu * Uydof
	eval_uz_gp = Nu * Uzdof

	#print '  // velocity(x) at gauss point //'
	#ccode_array(eval_ux_gp,1,1,'  u_gp')
	#print '  // velocity(y) at gauss point //'
	#ccode_array(eval_uy_gp,1,1,'  v_gp')
	#print '  // velocity(z) at gauss point //'
	#ccode_array(eval_uz_gp,1,1,'  w_gp')

	#ops = count_operations(eval_ux_gp)
	#ops = ops + count_operations(eval_uy_gp)
	#ops = ops + count_operations(eval_uz_gp)
	#total_ops = total_ops + ops

	dudx_gp = dNudx.transpose() * Uxdof
	dudy_gp = dNudy.transpose() * Uxdof
	dudz_gp = dNudz.transpose() * Uxdof

	dvdx_gp = dNudx.transpose() * Uydof
	dvdy_gp = dNudy.transpose() * Uydof
	dvdz_gp = dNudz.transpose() * Uydof

	dwdx_gp = dNudx.transpose() * Uzdof
	dwdy_gp = dNudy.transpose() * Uzdof
	dwdz_gp = dNudz.transpose() * Uzdof

		
	# compute pressure		
	Pdof = zeros((np,1))
	Np   = zeros((1,np))

	for i in range(0,np):
		Pdof[i,0] = X[nsd*nu+i]
		Np[0,i] = Symbol('Np['+str(i)+']')

	eval_p_gp  = Np * Pdof
	print '  // pressure at gauss point //'
	ccode_array(eval_p_gp,1,1,'  p_gp',False)

	ops = count_operations(eval_p_gp)
	total_ops = total_ops + ops
	
	# compute strain-rate
	eval_E_gp = zeros((ntens,1))
	eval_E_gp[0,0] = 0.5 * ( dudx_gp[0,0] + dudx_gp[0,0] )
	eval_E_gp[1,0] = 0.5 * ( dvdy_gp[0,0] + dvdy_gp[0,0] )
	eval_E_gp[2,0] = 0.5 * ( dwdz_gp[0,0] + dwdz_gp[0,0] )
	
	eval_E_gp[3,0] = ( dudy_gp[0,0] + dvdx_gp[0,0] )
	eval_E_gp[4,0] = ( dudz_gp[0,0] + dwdx_gp[0,0] )
	eval_E_gp[5,0] = ( dvdz_gp[0,0] + dwdy_gp[0,0] )


	print '  // strain-rate (e_xx,e_yy,e_zz,2.e_xy,2.e_xz,2.e_yz) at gauss point //'
	ccode_array(eval_E_gp,ntens,1,'  strain_rate_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_E_gp[d,0])

	total_ops = total_ops + ops

	E_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_E_gp[d,0] != 0:
			E_gp[d,0] = Symbol('strain_rate_gp['+str(d)+']')

	# compute constitutive law
	D = zeros((ntens,ntens))
	D[0,0] = 2.0 * eta_gp
	D[1,1] = 2.0 * eta_gp
	D[2,2] = 2.0 * eta_gp
	D[3,3] =       eta_gp
	D[4,4] =       eta_gp
	D[5,5] =       eta_gp

	# compute deviatoric stress
	eval_tau_gp = D * E_gp 
	print '  // deviatoric stress (t_xx,t_yy,t_zz,t_xy,t_xz,t_yz) at gauss point //'
	ccode_array(eval_tau_gp,ntens,1,'  tau_gp',False)
	ops = 0
	for d in range(0,ntens):
		ops = ops + count_operations(eval_tau_gp[d,0])

	total_ops = total_ops + ops

	tau_gp = zeros((ntens,1))
	for d in range(0,ntens):
		if eval_tau_gp[d,0] != 0:
			tau_gp[d,0] = Symbol('tau_gp['+str(d)+']')

	## Bt.D.B.u
	B = zeros((ntens,nsd*nu))
	for i in range(0,nu):
		B[0,nsd*i+0] = Symbol('dNudx['+str(i)+']')
		B[1,nsd*i+1] = Symbol('dNudy['+str(i)+']')
		B[2,nsd*i+2] = Symbol('dNudz['+str(i)+']')

		# exy = du/dy + dv/dx
		B[3,nsd*i+0] = Symbol('dNudy['+str(i)+']')
		B[3,nsd*i+1] = Symbol('dNudx['+str(i)+']')

		# exz = du/dz + dw/dx
		B[4,nsd*i+0] = Symbol('dNudz['+str(i)+']')
		B[4,nsd*i+2] = Symbol('dNudx['+str(i)+']')

		# eyz = dv/dz + dw/dy
		B[5,nsd*i+1] = Symbol('dNudz['+str(i)+']')
		B[5,nsd*i+2] = Symbol('dNudy['+str(i)+']')


	Bt_i_tau_gp = B.transpose() * tau_gp
	print '  // y = A11.u at gauss point //'
	#ccode_array(Bt_i_tau_gp,nsd*nu,1,'y0_uu')

	## -Bt.p.m
	m = zeros((ntens,1))
	m[0,0] = 1
	m[1,0] = 1
	m[2,0] = 1
	
	m[3,0] = 0
	m[4,0] = 0
	m[5,0] = 0

	p_gp = 0
	if eval_p_gp[0,0] != 0:
		p_gp = Symbol('p_gp')

	p_m_gp = p_gp * m

	grad_pressure_gp = -B.transpose() * p_m_gp
	print '  // y = A12.p at gauss point //'
	#ccode_array(grad_pressure_gp,nsd*nu,1,'y0_up')


	## Need to transpose this operator; -Bt.m.Np.p^h
	## -Np^T.m^T.B.u
	## -Np_i^T.m^T.(B.u)_j
	eval_div_gp = m.transpose() * E_gp
	print '  // divergence at gauss point //'
	ccode_array( eval_div_gp,1,1, '  div_gp', False )
	
	ops = count_operations(eval_div_gp)
	total_ops = total_ops + ops

	div_gp = 0
	if eval_div_gp[0,0] != 0:
		div_gp = Symbol('div_gp')

	Np_div_u = -Np.transpose() * div_gp
	
	print '  // y = A21.u at gauss point //'
	#ccode_array(Np_div_u,np,1,'y1_pu')


	print '  // reciprocal of viscosity average //'
	print '  one_on_eta_avg_gp = 1.0/avg_eta_gp;'
	total_ops = total_ops + 1

	print '  // shifted pressure at gauss point //'
	ones = zeros((np,1))
	for d in range(0,np):
		ones[d,0] = 1.0

	eval_sum_p = Pdof.transpose() * ones
	ccode_array( eval_sum_p,1,1,'  sum_p',False)
	ops = count_operations(eval_sum_p)
	total_ops = total_ops + ops

	sum_p = 0
	if eval_sum_p[0,0] != 0:
		sum_p = Symbol('sum_p')

	eval_scaled_p = -( (1.0/64.0) * sum_p )
	ops = count_operations(eval_scaled_p)
	total_ops = total_ops + ops

	#ccode_array( shifted_p_gp,1,1,'shifted_p_gp',False)
	print '  scaled_p = ' + str(eval_scaled_p) + ';'

	print '  scaled_p_on_eta_avg = scaled_p * one_on_eta_avg_gp;'
	print '  p_gp_on_eta_avg     = p_gp * one_on_eta_avg_gp;'
	total_ops = total_ops + 2


	scaled_p_on_eta_avg = Symbol('scaled_p_on_eta_avg')
	p_gp_on_eta_avg = Symbol('p_gp_on_eta_avg')

	stab = zeros((np,1))
	stab = - Np.transpose() * p_gp_on_eta_avg + ones * scaled_p_on_eta_avg 
#	ccode_array(stab,np,1,'  y1_pp',True)


	# combine results into one vector
	FAC = Symbol('FAC')
	for i in range(nu):
		Y[nsd*i  ] = FAC*(Bt_i_tau_gp[nsd*i  ,0] + grad_pressure_gp[nsd*i  ,0])
		Y[nsd*i+1] = FAC*(Bt_i_tau_gp[nsd*i+1,0] + grad_pressure_gp[nsd*i+1,0])
		Y[nsd*i+2] = FAC*(Bt_i_tau_gp[nsd*i+2,0] + grad_pressure_gp[nsd*i+2,0])

	for i in range(np):
#		Y[nsd*nu+i] = Np_div_u[i,0] + stab[i,0]
		Y[nsd*nu+i] = -FAC*( scaled_p_on_eta_avg + Np[0,i]*( div_gp + p_gp_on_eta_avg) )


	Yextract = extract_fields(Y,nsd,nu,vfields,np,pfield)
	lu = len(vfields)
	lp = len(pfield)
	
	ops = 0
	for d in range(0,lu*nu+lp*np):
		ops = ops + count_operations(Yextract[d,0])
	total_ops = total_ops + ops

	# add += 
	total_ops = total_ops + (lu*nu+lp*np)


	ccode_array(Yextract,lu*nu+lp*np,1,'  Y',True)
	print '  \n// total operations = ' + str(total_ops)

	print '}'
##########################################################################################		

def extract_fields(Y,nsd,nu,velfields,np,pressurefields):
#	print '// velfields = ',velfields
#	print '// pfields = ',pressurefields

	nvdofs = 0
	totalLength = 0
	if len(velfields) != 0:
		for f in velfields:
#			print '// fetching field...',f
			totalLength = totalLength + nu
			nvdofs = nvdofs + 1
	if len(pressurefields) != 0:
		totalLength = totalLength + np

#	print '// vdofs requested = ', nvdofs
#	print '// vector length = ', totalLength


	Yout = zeros((totalLength,1))
	cnt = 0
	if len(velfields) != 0:
		for i in range(0,nu):
			c = 0
			for f in velfields:
				if f == 'u':
					Yout[ nvdofs*i + c] = Y[ nsd*i + 0 ]
				if f == 'v':
					Yout[ nvdofs*i + c] = Y[ nsd*i + 1 ]
				if f == 'w':
					Yout[ nvdofs*i + c] = Y[ nsd*i + 2 ]
				c = c + 1

	if len(pressurefields) != 0:
		for i in range(0,np):
			Yout[nvdofs*nu + i] = Y[nsd*nu + i]


	return (Yout)


def zero_fields(X,nsd,nu,velfields,np,pressurefields):

	if len(velfields) != 0:
		for i in range(0,nu):
			for f in velfields:
				if f == 'u':
					X[ i + 0] = 0
				if f == 'v':
					X[ i + 1*nu] = 0
				if f == 'w':
					X[ i + 2*nu] = 0

	if len(pressurefields) != 0:
		for i in range(0,np):
			X[ i + nsd*nu ] = 0


def test3d():

	## BASIS FUNCTIONS
	nnodes_u = 27
	nnodes_p = 4
	dim = 3

	# velocity
	Udof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Udof[i,0] = Symbol('Ux['+str(i)+']')

	Vdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Vdof[i,0] = Symbol('Uy['+str(i)+']')

	Wdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Wdof[i,0] = Symbol('Uz['+str(i)+']')

	# pressure
	Pdof = zeros((nnodes_p,1))
	for i in range(0,nnodes_p):
		Pdof[i,0] = Symbol('P['+str(i)+']')

	print 'testing 3d'

	# zero components
#	for i in range(0,nnodes_u):
#		Udof[i,0] = 0
#		Vdof[i,0] = 0
#		Wdof[i,0] = 0
#	for i in range(0,nnodes_p):
#		Pdof[i,0] = 0

	X = zeros((nnodes_u*dim+nnodes_p,1))
	Y = zeros((nnodes_u*dim+nnodes_p,1))
	for i in range(0,nnodes_u):
		X[i         ,0] = Udof[i,0]
		X[i+nnodes_u,0] = Vdof[i,0]
		X[i+nnodes_u+nnodes_u,0] = Wdof[i,0]
	for i in range(0,nnodes_p):
		X[i+dim*nnodes_u,0] = Pdof[i,0]


	print '[input]'
	ccode_array(X,dim*nnodes_u+nnodes_p,1,'X',False)

	Stokes3dMixedFEM_MatMult_gp(nnodes_u,nnodes_p,X,Y)
	print '[[result y = A.x]]'
	ccode_array(Y,dim*nnodes_u+nnodes_p,1,'Y',False)

	print '[[result 2]]'
#	Yextract = extract_fields(Y,3,nnodes_u,['u','v','w'],nnodes_p,0)
#	ccode_array(Yextract,2*nnodes_u+0*nnodes_p,1,'Ye')

#	Yextract = extract_fields(Y,3,nnodes_u,0,nnodes_p,['p'])
#	ccode_array(Yextract,0*nnodes_u+1*nnodes_p,1,'Ye')
	
	Yextract = extract_fields(Y,3,nnodes_u,['u'],nnodes_p,['p'])
	ccode_array(Yextract,1*nnodes_u+1*nnodes_p,1,'Ye',True)

	
	print '\n\n*** FILTER y = A11.x ***'
	# zero components
	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp(nnodes_u,nnodes_p,Xc,Y)
	Yextract = extract_fields(Y,3,nnodes_u,['u','v','w'],nnodes_p,[])
	print '// Filtered result //'
	ccode_array(Yextract,3*nnodes_u+0*nnodes_p,1,'A11.u',True)


	print '\n\n*** FILTER y = A12.x ***'
	# zero components
	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v','w'],nnodes_p,[])
	Stokes3dMixedFEM_MatMult_gp(nnodes_u,nnodes_p,Xc,Y)
	Yextract = extract_fields(Y,3,nnodes_u,['u','v','w'],nnodes_p,[])
	print '// Filtered result //'
	ccode_array(Yextract,3*nnodes_u+0*nnodes_p,1,'A12.p',True)

	print '\n\n*** FILTER y = A21.x ***'
	# zero components
	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp(nnodes_u,nnodes_p,Xc,Y)
	Yextract = extract_fields(Y,3,nnodes_u,[],nnodes_p,['p'])
	print '// Filtered result //'
	ccode_array(Yextract,0*nnodes_u+1*nnodes_p,1,'A21.u',True)


	Xc = X * 1.0
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A',nnodes_u,nnodes_p,Xc,['u','v','w'],['p'])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A11',nnodes_u,nnodes_p,Xc,['u','v','w'],[])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v','w'],nnodes_p,[])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A12',nnodes_u,nnodes_p,Xc,['u','v','w'],[])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A21',nnodes_u,nnodes_p,Xc,[],['p'])


def test3d_2():

	## BASIS FUNCTIONS
	nnodes_u = 8
	nnodes_p = 8
	dim = 3

	# velocity
	Udof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Udof[i,0] = Symbol('Ux['+str(i)+']')

	Vdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Vdof[i,0] = Symbol('Uy['+str(i)+']')

	Wdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Wdof[i,0] = Symbol('Uz['+str(i)+']')

	# pressure
	Pdof = zeros((nnodes_p,1))
	for i in range(0,nnodes_p):
		Pdof[i,0] = Symbol('P['+str(i)+']')

	print 'testing 3d_2'

	X = zeros((nnodes_u*dim+nnodes_p,1))
	for i in range(0,nnodes_u):
		X[i         ,0] = Udof[i,0]
		X[i+nnodes_u,0] = Vdof[i,0]
		X[i+nnodes_u+nnodes_u,0] = Wdof[i,0]
	for i in range(0,nnodes_p):
		X[i+dim*nnodes_u,0] = Pdof[i,0]


	Xc = X * 1.0
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A',nnodes_u,nnodes_p,Xc,['u','v','w'],['p'])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A11',nnodes_u,nnodes_p,Xc,['u','v','w'],[])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v','w'],nnodes_p,[])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A12',nnodes_u,nnodes_p,Xc,['u','v','w'],[])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A21',nnodes_u,nnodes_p,Xc,[],['p'])


def generator_Q2Pm1_3d():

	## CLASS NAME
	MFOPClass = 'stokes_q2p1_mf_operators'

	## BASIS FUNCTIONS
	nnodes_u = 27
	nnodes_p = 4
	dim      = 3

	# velocity
	Udof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Udof[i,0] = Symbol('Ux['+str(i)+']')

	Vdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Vdof[i,0] = Symbol('Uy['+str(i)+']')

	Wdof = zeros((nnodes_u,1))
	for i in range(0,nnodes_u):
		Wdof[i,0] = Symbol('Uz['+str(i)+']')

	# pressure
	Pdof = zeros((nnodes_p,1))
	for i in range(0,nnodes_p):
		Pdof[i,0] = Symbol('P['+str(i)+']')
#	for i in range(0,nnodes_p):
#		Pdof[i,0] = 0



	X = zeros((nnodes_u*dim+nnodes_p,1))
	Y = zeros((nnodes_u*dim+nnodes_p,1))
	for i in range(0,nnodes_u):
		X[i         ,0] = Udof[i,0]
		X[i+nnodes_u,0] = Vdof[i,0]
		X[i+nnodes_u+nnodes_u,0] = Wdof[i,0]
	for i in range(0,nnodes_p):
		X[i+dim*nnodes_u,0] = Pdof[i,0]


	# write out the c file
	file = open(MFOPClass+'_def.c','w')
	sys.stdout = file

	print '\n\n\n#error(<<REMOVED AUTOGENERATED TAG>> ================== FILE ['+MFOPClass+'_def.c] ==================)\n\n'

	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_B',nnodes_u,nnodes_p,X,['u','v','w'],['p'])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_B11',nnodes_u,nnodes_p,Xc,['u','v','w'],[])


	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['v','w'],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_Buu',nnodes_u,nnodes_p,Xc,['u'],[])
	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','w'],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_Bvv',nnodes_u,nnodes_p,Xc,['v'],[])
	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v'],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_Bww',nnodes_u,nnodes_p,Xc,['w'],[])


	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v','w'],nnodes_p,[])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A12',nnodes_u,nnodes_p,Xc,['u','v','w'],[])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,[],nnodes_p,['p'])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A21',nnodes_u,nnodes_p,Xc,[],['p'])

	Xc = X * 1.0
	zero_fields(Xc,dim,nnodes_u,['u','v','w'],nnodes_p,[])
	Stokes3dMixedFEM_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_A22',nnodes_u,nnodes_p,Xc,[],['p'])

	file.close()
	sys.stdout = sys.__stdout__

def generator_Q2Pm1_3d_diag():

	## CLASS NAME
	MFOPClass = 'stokes_q2p1_mf_operators_diag'

	## BASIS FUNCTIONS
	nnodes_u = 27
	nnodes_p = 4
	dim      = 3

	# write out the c file
	file = open(MFOPClass+'_def.c','w')
	sys.stdout = file

	print '\n\n\n#error(<<REMOVED AUTOGENERATED TAG>> ================== FILE ['+MFOPClass+'_def.c] ==================)\n\n'

	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagB',nnodes_u,nnodes_p,['u','v','w'],['p'])


	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagB11',nnodes_u,nnodes_p,['u','v','w'],[])

	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagBuu',nnodes_u,nnodes_p,['u'],[])

	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagBvv',nnodes_u,nnodes_p,['v'],[])

	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagBww',nnodes_u,nnodes_p,['w'],[])

	Stokes3dMixedFEM_diag_MatMult_gp_with_filter('MatMultMF_Stokes_MixedFEM3d_diagA22',nnodes_u,nnodes_p,[],['p'])


	file.close()
	sys.stdout = sys.__stdout__


generator_Q2Pm1_3d()
generator_Q2Pm1_3d_diag()




