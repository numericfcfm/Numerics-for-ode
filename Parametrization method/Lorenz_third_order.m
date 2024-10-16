%Este programa se escribe sin acentos 
%Autor: Prof.Jaime Burgos Garcia
%Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
%Asignatura: Topicos Selectos de Analisis Numerico
%Tema: El metodo de parametrizacion para variedades invariantes de puntos de equilibrio.
% Caso del sistema de Lorenz
clc  
clear
format long
%valores propios
lambda1=-8/3;
lambda2=(1/2)*(-11-sqrt(1201));
%punto de equilibrio
p00=[0 0 0];
%vectores propios
p10=[0 0 1]; 
p01=[-(1/56)*(9+sqrt(1201)) 1 0];
P1= horzcat(p10',p01'); %termino de primer orden 
DF=[-10 10 0;28 -1 0;0 0 -8/3]; %Matriz Jacobiana
%%%%%%%%Orden 2
%%terminos de la derecha
p1320=P1(1,1)*P1(3,1);
p1311=P1(1,1)*P1(3,2)+P1(1,2)*P1(3,1);
p1302=P1(1,2)*P1(3,2);
p1220=P1(1,1)*P1(2,1);
p1211=P1(1,1)*P1(2,2)+P1(1,2)*P1(2,1);
p1202=P1(1,2)*P1(2,2);
%%%%solucion del sistema lineal para n=2 y m=0
A20=DF-(2*lambda1+0*lambda2)*eye(3);
s20=[0 p1320 -p1220];
p20=A20\s20';
%%%%solucion del sistema lineal para n=1 y m=1
A11=DF-(1*lambda1+1*lambda2)*eye(3);
s11=[0 p1311 -p1211];
p11=A11\s11';
%%%%solucion del sistema lineal para n=0 y m=2
A02=DF-(0*lambda1+2*lambda2)*eye(3);
s02=[0 p1302 -p1202];
p02=A02\s02';
%matriz de coeficientes de orden 2
P2= horzcat(p20,p11,p02)

%%%%%%%%%%%%%%%orden3
%%terminos de la derecha
p1330=P2(1,1)*P1(3,1)+P1(1,1)*P2(3,1);  
p1321=P2(1,1)*P1(3,2)+P2(1,2)*P1(3,1)+P1(1,1)*P2(3,2)+P1(1,2)*P2(3,1); 
p1312=P2(1,2)*P1(3,2)+P1(1,1)*P2(3,3)+P2(1,3)*P1(3,2)+P1(3,2)*P2(3,2); 
p1303=P2(1,3)*P1(3,2)+P1(1,2)*P2(3,3); 
p1230=P2(1,1)*P1(2,1)+P1(1,1)*P2(2,1); 
p1221=P2(1,1)*P1(2,2)+P2(1,2)*P1(2,1)+P1(1,1)*P2(2,2)+P1(1,2)*P2(2,1); 
p1212=P2(1,2)*P1(2,2)+P1(1,1)*P2(2,3)+P2(1,3)*P1(2,2)+P1(3,2)*P2(2,2); 
p1203=P2(1,3)*P1(2,2)+P1(1,2)*P2(2,3); 
%%%%solucion del sistema lineal para n=3 y m=0
A30=DF-(3*lambda1+0*lambda2)*eye(3);
s30=[0 p1330 -p1230]; 
p30=A30\s30';
%%%%solucion del sistema lineal para n=2 y m=1
A21=DF-(2*lambda1+1*lambda2)*eye(3);
s21=[0 p1321 -p1221]; 
p21=A21\s21';
%%%%solucion del sistema lineal para n=1 y m=2
A12=DF-(1*lambda1+2*lambda2)*eye(3);
s12=[0 p1312 -p1212]; 
p12=A12\s12';
%%%%solucion del sistema lineal para n=0 y m=3
A03=DF-(0*lambda1+3*lambda2)*eye(3);
s03=[0 p1303 -p1203]; 
p03=A03\s03';
%matriz de coeficientes de orden 3
P3= horzcat(p30,p21,p12,p03)

%codigoparaexportar de matlab a latex
%digits(14)
%latex_table = latex(vpa(sym(P2)))
%%%%%%%%%%Graficacion de la superficie
syms t1 t2
%parametrizacion a tercer orden
P=P1*[t1;t2]+P2*[t1^2;t1*t2;t2^2]+P3*[t1^3;t1^2*t2;t1*t2^2;t2^2];
%componentes
x=P(1);
y=P(2);
z=P(3);
%dominio de las variables locales
dominio=[-8.0 8.0 -8.0 8.0];
ezsurf(x,y,z,dominio)
title('Variedad de dimension dos')


