# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:00:27 2022

@author: Jaime Burgos García and Rey Alexis Salas Vega
"""
import numpy as np
import matplotlib.pyplot as plt
#import array_to_latex as a2l

#--- FUNCION QUE CALCULA CADA a_j -----
def A(p,j):
    s1=0
    for k in range(j+1):
        s1+=p[j-k,0]*p[k,0]
    return s1
#------------------------------------------
#--- FUNCION QUE CALCULA CADA s_{n-1} -----
def S(p,n): #A LA HORA DE LLAMAR PONER S(p,n-1)
    s2=0
    for i in range(1,n+1):
        s2+=p[0,0]*p[n+1-i,0]*p[i,0]+p[n+1-i,0]*A(p,i)
    return s2
#------------------------------------------
#------------PROGRAMA PRINCIPAL------------
pe=[0,0] #Coordenadas del punto de eq.
l=1 #Valor propio
#v=[1/np.sqrt(2),1/np.sqrt(2)] #vector propio normalizado
v=[1,1] #vector propio
o=20 #Orden del calculo
p=np.zeros([o+1,2]) #Genera la matriz con coeficientes de P(theta)
p[0:2,0:2]=np.array([pe,v])

for m in range(2,o+1):
    #print("paso: ",m)
    M=np.array([[-m*l,1],[1-3*p[0,0]**2,-m*l]])
    #print(M)
    b=np.array([0,S(p,m-1)])
    #print(b)
    z=np.linalg.solve(M,b)
    #print(z)
    p[m,:]=z
print(p)
#print(a2l.to_ltx(p,frmt='{:.8e}',arraytype='array'))
theta=np.linspace(0,2.5,101) #Parametro THETA
xx=np.zeros([1,101])
for i in range(101):
    s=0
    for j in range(o+1):
        s+=p[j,0]*theta[i]**j
    xx[0,i]=s
yy=np.zeros([1,101])
for i in range(101):
    s=0
    for j in range(o+1):
        s+=p[j,1]*theta[i]**j
    yy[0,i]=s
#print(xx)
#print(yy)
#plt.plot(xx[0,:],yy[0,:],'k')
#------------------------------------------
#--------------TRAYECTORIA EXACTA----------
t=np.linspace(-10,0,201)
x=np.sqrt(2)/np.cosh(t)
y=-np.sqrt(2)*np.sinh(t)/np.cosh(t)**2

#plt.axhline(0,color='k',lw=0.5)
plt.plot(x,y,'r',label='$\Gamma^+(t)$')
plt.xlim([0,1.5])
plt.ylim([0,0.8])
#plt.show()
#-----------------------------------------
plt.plot(xx[0,:],yy[0,:],'k--',label=r'Parameterization with $\theta \in[0,0.5]$')
plt.legend(loc="upper left",fontsize=9)
#plt.savefig('par-Duffing.pdf',bbox_inches='tight')
plt.show()
