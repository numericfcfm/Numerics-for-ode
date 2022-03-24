%Este programa se escribe sin acentos 
%Autor: Prof.Jaime Burgos Garcia
%Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
%Asignatura: Topicos Selectos de Analisis Numerico
%Tema: El metodo de parametrizacion para variedades invariantes de puntos de equilibrio.
% Caso de la ecuacion de Duffing
clear 
close all
clc
format long
%Bloque de parametros iniciales 
N=21;  %Numero que define el orden de la aproximacion deseada. El orden es N-1
p0=[0 0]; %coordenadas del punto de equilibrio
lambda=1;%eigenvalor para la parte inestable/estable
p1=[1 1]; %coordenadas del vector asociado a la parte inestable/estable
P=zeros(N,2); %arreglo para almacenar los coeficientes de la variedad
P(1,:)=p0;%coordenadas del punto de equilibrio en el primer renglon
P(2,:)=p1;%coordenadas del vector en el segundo renglon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 3:N %ciclo para calcular los coeficientes de orden 2 y superior
   

%%%%%%%%%%%%% calculo del termino s_n-1
sum1=0;
sum2=0;

 for j=2:n-1 %ciclo que calcula el segundo sumando de s_n-1
 sum3=0;
    for k=1:j %ciclo que calcula los terminos a_j
         sum3=sum3+P(n+1-k,1)*P(k,1);
         
     end
  sum1=sum1+P(n+1-j,1)*sum3;
end

for k=2:n-1%ciclo que calcula el primer sumando de s_n-1
    
  sum2=sum2+P(n+1-k,1)*P(k,1);
end
    s=P(1,1)*sum2+sum1; %termino s_n-1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%solucion de los sistemas lineales a cada paso
     M=[-(n-1)*lambda 1;1-3*(p0(1))^2 -(n-1)*lambda]; %matriz
     b=[0;s]; %vector derecho para el sistema
     P(n,:) = M\b; %solucion del sistema lineal
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

P %desplega el arreglo de coeficientes

%%%%%%%%bloque de graficacion de resultados 
xp=flipud(P(:,1))';
yp=flipud(P(:,2))'; %coeficientes en forma ascendente
thetamin=0;
thetamax=0.5; %extremos del dominio para theta
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
plot(xpv,ypv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%