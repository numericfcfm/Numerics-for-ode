%Este programa se escribe sin acentos 
%Autor: Prof.Jaime Burgos Garcia
%Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
%Asignatura: Topicos Selectos de Analisis Numerico
%Tema: El metodo de parametrizacion para variedades invariantes de puntos de equilibrio.
% Caso de la ecuacion de Duffing. VERSION CORREGIDA
%Advertencia. Este programa es una version corregida de la implementacion computacional del
%metodo presentado en el articulo, El metodo de parametrizacion para variedades invariantes de
%puntos de equilibrio de ecuaciones diferenciales ordinarias. Abstraction &
%Application, 30. pp 64-81, (2020). En consecuencia, los resultados
%proporcionados por este programa no corresponden a los presentados en
%dicho documento.
clear 
close all
clc
format long
%Bloque de parametros iniciales
%Por las necesidades de Matlab, el indice de las sumas comienzan desde 1.
N=21;  %Numero que define el orden de la aproximacion deseada. El orden es N-1
p0=[0 0]; %coordenadas del punto de equilibrio
lambda=1;%eigenvalor para la parte inestable
p1=[1 1]; %coordenadas d el vector asociado a la parte inestable
P=zeros(N,2); %arreglo para almacenar los coeficientes de la variedad
a=zeros(1,N); %arreglo para almacenar los coeficientes a_n (ver documento)
b=zeros(1,N); %arreglo para almacenar los coeficientes b_n (ver documento)
P(1,:)=p0;%coordenadas del punto de equilibrio 
P(2,:)=p1;%coordenadas del vector propio elegido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 3:N %ciclo para calcular los coeficientes de orden 2 y superior
   
sum1=0;
sum2=0;
%se ha usado el hecho que en el ejemplo actual el punto de equilibrio es el
%origen por lo cual los coeficientes a_1, a_2, b_1 y b_2 son cero (ver documento)

 for j=1:n %ciclo que calcula los elementos a_n
     
    sum1=sum1+P(n+1-j,1)*P(j,1);
           
 end
 a(n)=sum1;

for j=1:n %ciclo que calcula los elementos b_n
     
    sum2=sum2+P(n+1-j,1)*a(j);
           
end
 b(n)=sum2;
    s=sum2; %termino s_n-1, usando el hecho que p01=0 (ver documento)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%solucion de los sistemas lineales a cada paso
     M=[-(n-1)*lambda 1;1-3*(p0(1))^2 -(n-1)*lambda]; %matriz
     v=[0;s]; %vector derecho para el sistema
     P(n,:) = M\v; %solucion del sistema lineal
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

P %despliega el arreglo de coeficientes

%%%%%%%%bloque de graficacion de resultados 
xp=flipud(P(:,1))';
yp=flipud(P(:,2))'; %coeficientes en forma ascendente
thetamin=0;
thetamax=2.2; %extremos del dominio para theta
theta=thetamin:1/100:thetamax; %discretizacion del dominio
xpv=polyval(xp,theta);
ypv=polyval(yp,theta); %evaluacion de los polinomios dados por las coordenadas
tmin=0;
tmax=0.5; %extremos del dominio para t
t=-37:1/100:0; %discretizacion del dominio
x = sqrt(2)*sech(t);
y = -sqrt(2)*sech(t).*tanh(t);%evaluacion de la variedad en el dominio de t

figure(1)
plot(x,y)
hold on 
plot(xpv,ypv,'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%