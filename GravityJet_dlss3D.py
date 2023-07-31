#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 22:13:38 2021

@author: liuyiming
1 1 "Pbd_left"
1 2 "Pbd_right"
1 3 "vbd"
"""

"""
Dam Break by yiming liu

2021/5/28
"""

from dolfin import *
from dolfin.cpp.mesh import cells
from mshr import *
from ufl import indices
import matplotlib.pyplot as plt
#%matplotlib inline
import numpy as np
#from math import exp
import time

class CharacterLevelSetFunction(UserExpression):
	'''
	initialize the level set function
	Firstly, determine the position of an arbitrary point in the domain
	Secondly, calculate the distance to the minimum distance to the interface
	Finally, calculate the CharacterLevelSetFunction, 
	in the interface =1, out the interface = 0, on the transition area of interface it is continuous
	'''
	def __init__(self, param1, param2, param3, param4, **kwargs):
		super().__init__(**kwargs)
		self.param1 = param1 #radius
		self.param2 = param2 #middle point x coordinate value
		self.param3 = param3 #epsilon
		self.param4 = param4 #middle point y coordinate value
	'''
	def eval(self, values, x):
		dist=pow( pow(x[0]-self.param2, 2) + pow(self.param4- x[1], 2), 0.5)
		if dist <= self.param1:# and x[0]<=0.5:
			phi_d = dist - self.param1
			values[0] = 1/(1+exp(phi_d/self.param3))
		if dist > self.param1:# and x[0]<=0.5:
			phi_d = dist - self.param1
			values[0] = 1/(1+exp(phi_d/self.param3))
	'''
	def eval(self, values, x):
		R_air=pow( pow(x[0], 2) + pow(x[2], 2), 0.5)
		dist_out = pow( pow(R_air - self.param1,2) + pow(self.param4-1-x[1],2),0.5)
		if x[1] > 10.:# and x[0]<=0.5:
			phi_d = 10. - x[1]
			values[0] = 1/(1+exp(phi_d/self.param3))
		elif R_air <= self.param1:
			phi_d = 10.  - x[1]
			values[0] = 1/(1+exp(phi_d/self.param3))
		elif R_air > self.param1:
			phi_d = dist_out
			values[0] = 1/(1+exp(phi_d/self.param3))
processID =  MPI.comm_world.Get_rank()
set_log_level(40) #50, 40, 30, 20, 16, 13, 10


markvalue2= 0 #1-variable viscosity, 0-constant vis


res=12.5 #resolution: set the number following by the fname which depends on the selected mesh

#res=10 #resolution
scale=1
t_end = 10#8
B=0.0125 #0.001*res
dt = B*(scale/res)
t=0.
left_x=-4
right_x=4
bottom=0
top=10
pip_d=1
pip_h=1


foldername = "data/example1/"


fname = "geofolder/Jetflow_3D_res12_5"#resolution is 12.5

#fname = "geo/Jetflow_3D_uni"
mesh = Mesh(fname + '.xml')

open_d=pip_d
bd_left = left_x
bd_right = right_x
pressurebd = top
r=0.5*pip_d
xc=0
yc=top+pip_h
#%%
i, j, k, l = indices(4)
delta = Identity(3)

scalar = FiniteElement('P', mesh.ufl_cell(), 1) #pressure 
vector = VectorElement('P', mesh.ufl_cell(), 1) #velocity
vector1 = VectorElement('P', mesh.ufl_cell(), 1) #interface normal
scalar1 = FiniteElement('P', mesh.ufl_cell(), 1) #level set function
mixed_element = MixedElement([scalar, vector])
Space = FunctionSpace(mesh, mixed_element)
ScalarSpace = FunctionSpace(mesh, scalar) #pressure
ScalarSpace1 = FunctionSpace(mesh, scalar1) #level set function
VectorSpace = FunctionSpace(mesh, vector) #velocity
VectorSpace1 = FunctionSpace(mesh, vector1) #interface normal

## Define boundary condition
no_slip = "on_boundary"
pbd = "on_boundary && near(x[1], pressurebd)  "
piptop = "on_boundary && near(x[1], yc) "
bdbottom = "on_boundary && near(x[1], bottom) "
outlet = "on_boundary && near(x[1], 0) && x[0]<=0.25 && x[0]>= -0.25 "
sub_domains = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1)
sub_domains.set_all(0)
noslip = CompiledSubDomain(no_slip)
noslip.mark(sub_domains, 1)
pbd = CompiledSubDomain(pbd, pressurebd=pressurebd)
pbd.mark(sub_domains,2)
piptop = CompiledSubDomain(piptop, yc=yc)
piptop.mark(sub_domains,3)
bdbottom = CompiledSubDomain(bdbottom, bottom=bottom)
bdbottom.mark(sub_domains,4)
di = Measure('dS', domain=mesh, subdomain_data=sub_domains)
da = Measure('ds', domain=mesh, subdomain_data=sub_domains)
dv = Measure('dx', domain=mesh, subdomain_data=sub_domains)
#cutopen = CompiledSubDomain(cut_open, y_c=cut_off,xmid=middle_point_x,d=open_d)
#cutopen.mark(sub_domains, 2)

#file = File(foldername + "subdomains.pvd")
#file << sub_domains
# set parameters : x_ref = 10mm g_ref = 9.8E3mm/s2 u_ref = (x_ref*g_ref)^0.5 = 313 mm/s
#rho_ref =10E-3 #g/mm3 rho_liquid = 10E-3 #g/mm3  rho_gas = 10 E-6.: rho0 =1 rho1 =0.001
# mu1_real = 0.1 Pa s/ 0.001 Pa s
#shear_rate_ref = u_ref/x_ref = 10/313
n=0.599 
g_ref = 9.8E3   #mm/s2
sig_ref = 72.7 #72.7    #muN/mm
rho_ref =1.26E-3  #g/mm3
x_ref = 8.    #mm
u_ref = 200#(x_ref*g_ref)**0.5 #mm/s
shear_rate_ref = (u_ref/x_ref)**(n-1) #/s
K_ref=1000.        #Pa s2
mu_ref= 1.3 #K_ref*shear_rate_ref
p_ref = rho_ref*x_ref*g_ref
#%%
bc = []
bcv = []
bcp = []
bcn = []
bcphi = []
# No-slip condition
Amp = 1.#70./u_ref
#Amp = u_ref
p_initial = 0#1E5/p_ref
#inflow_profile = Expression( ('0', 'Amp*(x[0]-bd_left)*(x[0] + bd_left)'), degree=2, Amp=Amp,bd_left=bd_left)
bc.append(DirichletBC(Space.sub(1), Constant((0.,0.,0.)), noslip))
bc.append(DirichletBC(Space.sub(1), Constant((0.,-Amp,0.)), piptop))
bc.append(DirichletBC(Space.sub(0), p_initial, pbd))
bcphi.append(DirichletBC(ScalarSpace1, Constant(1.), piptop))


dunkn = TrialFunction(Space)
test = TestFunction(Space)
unkn = Function(Space)
unkn0 = Function(Space)
unkn00 = Function(Space)
p0, u0 = split(unkn0)
p, u = split(unkn)
q, v = split(test)
n_phi = Function(VectorSpace1)
del_n_phi = TestFunction(VectorSpace1) 

del_phi_c = TestFunction(ScalarSpace1) #level set function
phi_c = Function(ScalarSpace1) #level set function
phi_cr = Function(ScalarSpace1) #level set function after reinitialization
phi_c0 = Function(ScalarSpace1)
phi_cr0 = Function(ScalarSpace1) #level set function after reinitialization


V = FunctionSpace(mesh, 'P', 1)


## Define initial conditions
C=1.# narrow_100-1.5, wide_100-1.1, narrow_200-3, wide_200-2.1 regular_narrow_100-1.3 regular_wide_100_1.1 curve_irregular_200:2.5 

epsilon=C*pow(scale/res,0.9)#related to the size of element
#dts=C*pow(scale/res,1.1) #sub-time step related to the mesh resolution
dts=dt #sub-time step related to the mesh resolution

initial_phi_c = CharacterLevelSetFunction(param1=r, param2=xc, param3=epsilon, param4=yc, element=V.ufl_element(),degree=1)#initial level set function
initial_v = Expression(('0.', '0.', '0.'), degree=1)

#%%
## initialization
#vv0 = project(initial_v, VectorSpace)
phi_c.assign(project(initial_phi_c, ScalarSpace1))
phi_c0.assign(project(initial_phi_c, ScalarSpace1))
#assign(unkn.sub(1), vv0)
#assign(unkn0.sub(1), vv0)
#mutiply the inital level-set value can make the initial calculation step more stable, but it has no effect on the result
# can also use p0.assign(project(initial_p, ScalarSpace)) 



Re=rho_ref*u_ref*x_ref/mu_ref
We=rho_ref*u_ref*u_ref*x_ref/sig_ref
Fr=u_ref/(x_ref*g_ref)**0.5

if processID == 0: print('Re= %.8s'%Re,'We= %.6s'%We,'Fr= %.6s'%Fr,'mu= %.6s'%mu_ref)
K0=1.#K_ref#1.
#K0=K_ref#1.
mu0=1.#mu_ref#1.
#mu0=mu_ref#1.
rho0=1.
rho1=1E-3
#rho0= rho_ref
#rho1=1E-6
mu1=1E-3/mu_ref #1E-3
#mu1=1E-3 #1E-3
# write file
#file_pp = File(foldername+'/p_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.pvd'%Fr)
#file_vv = File(foldername+'/v_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.pvd'%Fr)
#file_c = File(foldername+'/c_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.pvd'%Fr)
#file_n = File(foldername+'/n_phi_solve.pvd')
#file_n = File(foldername+'/n_phi_project.pvd')
#file_c << (phi_c0,t)
#file_pp << (p0,t)

file_p = XDMFFile(foldername+'/p_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.xdmf'%Fr)
file_v = XDMFFile(foldername+'/v_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.xdmf'%Fr)
file_phi = XDMFFile(foldername+'/phi_Re%.6s' %Re + '_dt%.5s'%dt + '_res%d'%res + '_We%.3s'%We+ '_Fr%.3s.xdmf'%Fr)
file_p.write_checkpoint(unkn.sub(0), "p", t, append=False)
file_v.write_checkpoint(unkn.sub(1), "v", t, append=False)
file_phi.write_checkpoint(phi_c, "phi", t, append=False)


g_= 1.
#g_= g_ref
rou= phi_c*rho0+(1-phi_c)*rho1
rou_0= phi_c0*rho0+(1-phi_c0)*rho1

ndim = mesh.geometry().dim()
n_bd = FacetNormal(mesh)
# viscous operator
def a(phi,chi, psi): return Constant(0.5) * inner(grad(phi) + grad(chi).T, grad(psi) + grad(psi).T)
# divergence operator
def b(phi, psi): return inner(div(phi), psi)
# non-linear convection operator
def c(phi, chi, psi): return dot(dot(grad(chi), phi), psi)
# Crank-Nicholson schema average
def d(phi, psi): return Constant(0.5) * (phi + psi)
def e(phi, psi): return dot(grad(phi), psi)
# interface normal n = ▽ϕ/|ϕ|
#normalphi=grad(phi_c)/pow(dot( grad(phi_c), grad(phi_c) ), 0.5)
nphi = grad(phi_c)/pow(dot( grad(phi_c), grad(phi_c) ), 0.5)
F_n_phi = (dot(grad(phi_c)/pow(dot( grad(phi_c), grad(phi_c) ), 0.5), del_n_phi) - dot(n_phi,del_n_phi)) * dx
#n_phi.assign(project(nphi,VectorSpace1, 'gmres'))
g = Expression(('0','-g','0.'),degree=1,g=g_)
T = ( Identity(len(n_phi)) - outer(n_phi, n_phi) ) * pow(dot( grad(phi_c), grad(phi_c) ), 0.5)  #suface tension fs = ▽·T
#T = ( Identity(len(nphi)) - outer(nphi, nphi) ) * pow(dot( grad(phi_c), grad(phi_c) ), 0.5)  #suface tension fs = ▽·T

d_ = as_tensor(1.0/2.0*(u[k].dx(l)+u[l].dx(k)) , (k,l))
d_dev = as_tensor(d_[i,j] - 1./3.*d_[k,k]*delta[i,j] , (i,j))
II = as_tensor(d_dev[i,j]*d_dev[i,j]+0.00001, ())

#
if markvalue2 == 1: mu= phi_c*K0*II**(0.5*(n-1)) + (1-phi_c)*mu1
if markvalue2 == 0: mu= phi_c*mu0+(1-phi_c)*mu1
sig = Constant(sig_ref)
tau = as_tensor(  mu*d_dev[i,j] , (i,j) )
Is = Identity(len(nphi)) - outer(nphi, nphi)
dlta = pow(dot( grad(phi_c), grad(phi_c) ), 0.5)

#  momentum balance eq
F_momentum = (rou*dot((u-u0),v)/dt + rou*c(u,u,v) + (1. / Re) * inner(tau,grad(v)) + e(p, v) - (1./(Fr*Fr))*rou*dot(g,v) + (1./We)*inner(T,grad(v)) ) * dv 
F_mass = b(u,q)*dv
F_affixation = (dot((u-u0),grad(q)) - (dt/(Fr*Fr))*dot(g,grad(q)) + dt*c(u,u,grad(q)) + (dt/rou)*e(p, grad(q)) - (dt/rou)*(1./Re)*dot(div(tau),grad(q))- (dt/rou)*(1./We)*dot(div(T),grad(q)) ) * dv

F=F_momentum+F_mass+F_affixation
#F=F_momentum+F_mass

## advection and reinitialization Crank-Nicholson schema
F_advection = ( (1/dt)*dot(phi_c-phi_c0, del_phi_c) - 0.5*dot( (phi_c+phi_c0), dot( u, grad(del_phi_c) ) ) ) * dv #+ 0.5*(phi_c+phi_c0)*dot(u,n_bd)*del_phi_c*(da(2)+da(4))

F_reinitialization = ( (1/dts)*dot(phi_cr - phi_cr0, del_phi_c) - dot(d(phi_cr, phi_cr0)*(1-d(phi_cr, phi_cr0)), dot(n_phi, grad(del_phi_c))) + epsilon*dot( grad(d(phi_cr, phi_cr0)), grad(del_phi_c))  ) * dv #+ dot(d(phi_cr, phi_cr0)*(1-d(phi_cr, phi_cr0))*n_phi - epsilon*grad(d(phi_cr, phi_cr0)), n_bd)*del_phi_c*(da(2)+da(4))

t0=0
t1=0
t2=0
t3=0
Gain1 = derivative(F, unkn, dunkn)
problem_NS = NonlinearVariationalProblem(F, unkn, bc, Gain1)
solver_NS  = NonlinearVariationalSolver(problem_NS)
tol1=1E-3
tol2=1E-5
tol3=1E-4
tol4=1E-6
solvetype = "gmres" #gmres
solver_NS.parameters["newton_solver"]["relative_tolerance"] = tol1
solver_NS.parameters["newton_solver"]["absolute_tolerance"] = tol2
solver_NS.parameters["newton_solver"]["convergence_criterion"] = "residual"
solver_NS.parameters["newton_solver"]["error_on_nonconvergence"] = True
solver_NS.parameters["newton_solver"]["linear_solver"] = solvetype
solver_NS.parameters["newton_solver"]["lu_solver"]["symmetric"] = False 
solver_NS.parameters["newton_solver"]["maximum_iterations"] = 1000
#solver_NS.parameters["newton_solver"]["preconditioner"] = "amg"
solver_NS.parameters["newton_solver"]["relaxation_parameter"] = 1.0
solver_NS.parameters["newton_solver"]["report"] = False

solver_NS.parameters["newton_solver"]["krylov_solver"]["nonzero_initial_guess"] = True #starting from the sloution at the last time step
solver_NS.parameters["newton_solver"]["krylov_solver"]["relative_tolerance"] = tol3
solver_NS.parameters["newton_solver"]["krylov_solver"]["absolute_tolerance"] = tol4
solver_NS.parameters["newton_solver"]["krylov_solver"]["monitor_convergence"] = False
solver_NS.parameters["newton_solver"]["krylov_solver"]["maximum_iterations"] = 10000

krylov_params = parameters["krylov_solver"]
krylov_params["nonzero_initial_guess"] = True
krylov_params["relative_tolerance"] = tol3
krylov_params["absolute_tolerance"] = tol4
krylov_params["monitor_convergence"] = False
krylov_params["maximum_iterations"] = 10000
#%%
nn=0

start_time=time.time()
while t < t_end:
	t += dt
	nn+= 1
	dofs=len(unkn.vector())
	if processID == 0: print('time: %f s,  with %d dofs' % (t,dofs) )
	
	begin("Computing Level Set Fucntion")
	LS_start=time.time()
	try:
		solve(F_advection==0, phi_c, bcphi, solver_parameters={"newton_solver":{"linear_solver": solvetype, "convergence_criterion": "residual","relative_tolerance": tol1,"absolute_tolerance":tol2,\
			"krylov_solver": krylov_params, "maximum_iterations":1000, "error_on_nonconvergence": True} })
	except:
		if processID == 0: print('Advaction not converged --------------------------------------' )
		solve(F_advection == 0, phi_c, bcphi, solver_parameters={"newton_solver":{"linear_solver": solvetype, "preconditioner": "none", "relative_tolerance": tol1, "relaxation_parameter":0.9, "maximum_iterations":1000, "error_on_nonconvergence": False} } )
	LS_cost = time.time() - LS_start
	t2 += LS_cost
	end()

	nreinitial=0
	maxstep=1000
	tol=3*epsilon/C
	#phi_cr0.assign(phi_c)
	assign(phi_cr0, phi_c)
	#n_phi.assign(project(nphi,VectorSpace1))
	solve(F_n_phi==0, n_phi, bcn, solver_parameters={"newton_solver":{"linear_solver": solvetype, "convergence_criterion": "residual","relative_tolerance": tol1,"absolute_tolerance":tol2,\
		"krylov_solver": krylov_params, "maximum_iterations":1000, "error_on_nonconvergence": True} })
	begin("Reinitialization")
	Re_start = time.time()
	while nreinitial < maxstep:
		nreinitial += 1
		solve(F_reinitialization==0, phi_cr, bcphi, solver_parameters={"newton_solver":{"linear_solver": solvetype, "convergence_criterion": "residual","relative_tolerance": tol1,"absolute_tolerance":tol2,\
			"krylov_solver": krylov_params, "maximum_iterations":1000, "error_on_nonconvergence": True} })
		#error=(phi_cr.vector() - phi_cr0.vector()).max() / dts
		error = assemble(abs(phi_cr-phi_cr0)*dv)/dts
		#error = norm(phi_cr - phi_cr0)/dts
		#error = abs((phi_cr.vector() - phi_cr0.vector())[0:]).max() / dts
		if error < tol and nreinitial > 1: break
		if processID == 0:
			if nreinitial==1: print('Error at substep 1: %7f' % error)
			
		#phi_cr0.assign(phi_cr)
		assign(phi_cr0, phi_cr)

	Re_cost = time.time()-Re_start
	t3 += Re_cost
	end()
	if processID == 0: print('Error in the end: %7f' % error)
	if processID == 0: print('Number of iterations for reinitialization: %d' % nreinitial)

	#phi_c.assign(phi_cr)
	#unkn0.assign(unkn)
	#phi_c0.assign(phi_c)
	assign(phi_c,phi_cr)
	#n_phi.assign(project(nphi,VectorSpace1))
	'''
	solve(F_n_phi==0, n_phi, bcn, solver_parameters={"newton_solver":{"linear_solver": solvetype, "convergence_criterion": "residual","relative_tolerance": tol1,"absolute_tolerance":tol2,\
		"krylov_solver": krylov_params, "maximum_iterations":1000, "error_on_nonconvergence": True} })
	'''
	begin("Computing N-S")
	NS_start=time.time()
	try:
		solver_NS.solve()
	except:
		if processID == 0: print('NS not converged --------------------------------------' ) 
		solve(F == 0, unkn, bcs=bc, J=Gain1, solver_parameters={"newton_solver":{"linear_solver": solvetype, "preconditioner": "none", "relative_tolerance": tol1, "relaxation_parameter":0.9, "maximum_iterations":1000, "error_on_nonconvergence": False} } )
	#p, u = unkn.split(deepcopy=True)
	NS_cost=time.time() - NS_start
	t1 += NS_cost
	end()
	assign(unkn0,unkn)
	assign(phi_c0,phi_c)
	if nn == 1 or (nn % 200 == 0 and nn<=10000000):
		div_V=assemble(div(rou*u)*dv)
		if processID == 0: print('div_u= %s' %div_V) 
		file_p.write_checkpoint(unkn.sub(0), "p", t, append=True)
		file_v.write_checkpoint(unkn.sub(1), "v", t, append=True)
		file_phi.write_checkpoint(phi_c, "phi", t, append=True)
		#file_pp << (p0,t)
		#file_c << (phi_c0,t)
		#file_vv << (u_c,t)
		#file_n << (n_phi,t)
total_time = time.time()-start_time
if processID == 0: print('Re= %.6s'%Re,'We= %.2s'%We,'mu= %.6s'%mu_ref)
if processID == 0: print("Projection process taken: %.3f s" %t1)
if processID == 0: print("Solve level set function taken: %.3f s" %t2)
if processID == 0: print("Reinitialization process taken: %.3f s" %t3)
if processID == 0: print("total time taken: %.3f s" %total_time)
div_V=assemble(div(rou*u)*dv)
if processID == 0: print('div_u= %s' %div_V) 
#file_pp << (p0,t)
#file_c << (phi_c0,t)
#file_vv << (u_c,t)
