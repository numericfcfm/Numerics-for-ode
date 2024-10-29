# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 08:42:10 2024
Este programa se escribe sin acentos 
@author: Rey Alexis Salas Vega and Jaime Burgos Garcia

Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
Asignatura: Topicos Selectos de Analisis Numerico
Tema: El metodo de parametrizacion para la variedad estable de dos dimensiones del origen
en el sistema de Lorenz 
Antes de ejecutar el archivo lea el articulo, El metodo de parametrizacion para variedades invariantes de
puntos de equilibrio de ecuaciones diferenciales ordinarias. Abstraction &
Application, 30. pp 64-81, (2020)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
#import array_to_latex as a2l
import sympy as sp
from matplotlib import cm


# Vector Field of Lorenz system
sigma=10
rho=28
beta=8/3

def f(t,Y):
    x,y,z=Y
    return sigma*(y-x), rho*x-x*z-y, x*y-beta*z

# Subprogram S_{mn}
def S(p1,p2,m,n):
    return np.sum(p[p1,:m+1,:n+1]*np.flip(p[p2,:m+1,:n+1]))


ep=np.array([0,0,0]) # Equilibrium point
l1=-8/3 # Eigenvalues
l2=(-11-np.sqrt(1201))/2
ev1=np.array([0,0,1]) # Stable eigenvectors
ev2=np.array([(-9-np.sqrt(1201))/56,1,0])

############### Integration of the unstable manifold
epsilon=1e-1
ev3=np.array([(-9+np.sqrt(1201))/56,1,0]) #unstable eigenvectors
iv=ep+ev3*epsilon #initial value
int=np.array([0,50]) #interval of t
t=np.linspace(0,50,10001) #specific points of the solution

sol=solve_ivp(f,int, iv, method='DOP853',t_eval=t,rtol=1e-13,atol=1e-14) #DOP853
XX=sol.y[0,:]
YY=sol.y[1,:]
ZZ=sol.y[2,:]

############### Plotting the one-dimensional invariant manifold
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.view_init(elev=-2, azim=-75)#, roll=15)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')

ax.plot(XX,YY,ZZ,color='r',label='Unstable manifold',linewidth=0.3)

################## Parameterization method 
o=20 # Order
p=np.zeros([3,o+1,o+1]) #Coefficient tensor
p[:,0,0]=ep #Terms of order 0
p[:,1,0]=ev1 #Terms of order 1
p[:,0,1]=ev2
#print(p)

# Construction of the symbolic polynomial
t1=sp.Symbol('t1') #variable theta1 as in page 74 of the above-referenced work
t2=sp.Symbol('t2') #variable theta1 as in page 74 of the above-referenced work
par=sp.Array(ep+ev1*t1+ev2*t2)

#Loops for the computation of higher order terms
for i in range(2,o+1):
    for m in range(0,i+1):
        M=np.array([
            [-sigma, sigma, 0],
            [rho, -1, 0],
            [0, 0, -beta]
            ])-(m*l1+(i-m)*l2)*np.eye(3)
        b=np.array([0,S(0,2,m,i-m),-S(0,1,m,i-m)])
        a=np.linalg.solve(M,b)
        p[:,m,i-m]=a
        # Add terms to Parameterization
        par+=sp.Array(a*t1**m*t2**(i-m))

#print(p)
#print(par)


# Domain to plot the surface
pts=[81,8] # Dimension mesh 
# the variables t1 and t2 in polar coordinates (r,alpha)
def P(alpha,r):
    return np.array(par.subs(t1,r*np.cos(alpha)).subs(t2,r*np.sin(alpha)))

theta1=np.zeros(pts[0])
rmax=15 #Size of the fundamental domain
thetamax=2*np.pi
theta1=np.linspace(0,thetamax,pts[0]) # alpha
theta2=np.linspace(0,rmax,pts[1]) # r

T1,T2=np.meshgrid(theta1,theta2)

xx,yy,zz=np.zeros([3,pts[1],pts[0]])

#loop to obtain the evaluation of the parametrization on each point of the mesh
for i in range(0,pts[1]):
    for j in range(0,pts[0]):
        print('[',i,',',j,']')
        print(P(T1[i,j],T2[i,j]))
        xx[i,j],yy[i,j],zz[i,j]=P(T1[i,j],T2[i,j])

## Plotting the two-dimensional invariant manifold

surf=ax.plot_surface(xx,yy,zz,cmap=cm.autumn,alpha=1)#Other styles: seismic,viridis,plasma,hot,gist_rainbow,winter,cool,copper,binary
ax.set_axis_off()
ax.set_facecolor('black')
plt.show()
