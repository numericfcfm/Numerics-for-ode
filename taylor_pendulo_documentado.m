%Metodo de Taylor con diferenciacion automatica para problemas de valor
%inicial con paso constante. Caso del pendulo simple
%Este programa se escribe sin acentos 
%Autor: Jaime Burgos Garcia
%Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
clear
clc
format long
%%%%bloque de condiciones iniciales
u10=3.14; %anguloinicial
u20=0;    %velocidad inicial
u30=sin(u10); %condicion inicial para variable auxiliar
u40=cos(u10); %condicion inicial para variable auxiliar
%%%%%%%%%%%%%%%%%%%%%%%
p=13;%p definira el orden del integrador que sera p-1
k=sin(u10/2);
m=k^2;
T=4*ellipticK(m); %periodo dado por una funcion eliptica con elementos k y m
tmin=0;
tmax=T;
N=160; %numero de nodos requeridos en la integracion
h=(tmax-tmin)/N; %tamano de paso constante
w=zeros(N,4); % arreglo de Nx4 para guardar las aproximaciones en cada nodo
t=zeros(1,N);  % arreglo de Nx1 para guardar la discretizacion en el tiempo
t(1)=tmin;
w(1,:)=[u10;u20;u30;u40]' ;%condiciones iniciales
tic %instruccion tic-toc para medir el tiempo de ejecucion
for i = 1:N
%%%%%%%%%%calculo de los coeficientes (derivadas normalizadas) en cada nodo t_i  
u1=zeros(1,p); %defincion de 4 arreglos que contendran los coeficientes de la serie de cada coordenada
u1(1)=w(i,1);
u2=zeros(1,p);
u2(1)=w(i,2);
u3=zeros(1,p);
u3(1)=w(i,3);
u4=zeros(1,p);
u4(1)=w(i,4);
 %%%%%%%%%%%%%%% ciclos para calcular los coeficientes por medio de
 %%%%%%%%%%%%%%% formulas recursivas
     for n=1:p-1
    u1(n+1)=u2(n)/(n);
    u2(n+1)=-u3(n)/(n);
    
    sum=0;
    for r=1:n
     sum=sum+u2(n+1-r)*u4(r);
     end
    u3(n+1)=sum/(n);
    
    sum=0;
    for r=1:n
     sum=sum+u2(n+1-r)*u3(r);
     end
    u4(n+1)=-sum/(n);
   
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%fin del calculo de los coeficientes
   %%%%%%%%%%%%%%%%%%%%%%% ciclos para calcular las aproximaciones w_i
    w1=0;
    for j=1:p
       w1=w1+u1(j)*h^(j-1);
    end
    w(i+1,1)=w1;
      
      w2=0;
    for j=1:p
       w2=w2+u2(j)*h^(j-1);
    end
      w(i+1,2)=w2;
      
      w3=0;
    for j=1:p
       w3=w3+u3(j)*h^(j-1); 
    end
      w(i+1,3)=w3;
      
      w4=0;
    for j=1:p
       w4=w4+u4(j)*h^(j-1); 
    end
       w(i+1,4)=w4;
      t(i+1)=i*h; %actualizacion del tiempo
  
end
%%%%%fin de la rutina
toc
figure(1) 
plot(w(:,1),w(:,2)) %graficacion de la orbita aproximada obtenida
tau=ellipticK(m)-t';
[s,c,d] = ellipj(tau,m); %funcion especial de Jacobi
theta=2*asin(k*s); %evaluacion de la solucion analitica en cada nodo
%abs(w(:,1)-theta)% %error absoluto en cada paso (el angulo theta se vuelve cero)
abs(w(end,1)-theta(end))/abs(theta(end)) % error relativo en el punto final