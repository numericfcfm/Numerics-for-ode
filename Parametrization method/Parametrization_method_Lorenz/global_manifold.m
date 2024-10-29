%Este programa se escribe sin acentos 
%Autor: Prof.Jaime Burgos Garcia
%Facultad de Ciencias Fisico Matematicas. Universidad Autonoma de Coahuila
%Asignatura: Topicos Selectos de Analisis Numerico
%Tema: Esquema basico para la extension de la variedad estable local de dos dimensiones del origen
%en el sistema de Lorenz 
%Antes de ejecutar el archivo lea el articulo, El metodo de parametrizacion para variedades invariantes de
%puntos de equilibrio de ecuaciones diferenciales ordinarias. Abstraction &
%Application, 30. pp 64-81, (2020)
%SE UTILIZO UNA VARIEDAD ESTABLE LOCAL DE ORDEN 20 CON DOMINIO FUNDAMENTAL
%DE RADIO 16
clear 
clc
%%%%%%%%%%%%Datos iniciales
ic=load('border.m'); %se cargan los datos provenientes de la discretizacion del borde de la variedad estable 
Tmax=150; %Distancia de longitud de arco
dim=size(ic);
N=dim(1);
tspan=[Tmax 0]; %intervalo de integracion hacia atras
options = odeset('RelTol',2.22045e-014,'AbsTol',eps); %opcion de precision de la integracion

%%%%%%%%%%Integracion numerica de cada punto en el borde
for k=1: N
vec=ic(k,:);%condicion inicial k-esima
[t1,L]=ode113(@lorenzfield_normalized,tspan,vec,options);
figure(1)
plot3(L(:,1),L(:,2),L(:,3),'Color','b','LineWidth',0.5)
hold on
%view(60,30) %vista de la solucion graficada
xlabel('x')
ylabel('y')
zlabel('z')
end
%Graficacion de todas las trayectorias
hold off