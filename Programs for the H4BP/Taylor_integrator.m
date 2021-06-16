%Numerical integrator for the Hill four-body problem
%Author: Jaime Burgos-Garcia. FCFM. Autonomous University of Coahuila
clear 
clc   
format long
%%%%initial conditions
m=0.5;  %value of the mas parameter \mu
%initial positions
x0=1.022159719609811 ;
y0=0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
z0=0;
%initial velocities
xdot0=0; 
ydot0=-1.392074692111009;  
zdot0=0.587678513478811 ;   
period= 2*5.090815423831103; %period of the orbit                       
tmin=0;       
tmax=period; % Integration time;
x70=1/sqrt(x0^2+y0^2+z0^2); %computing the auxiliar initial condition
x=[x0 y0 z0 xdot0 ydot0 zdot0 x70]; %seven dimensional vector of initial conditions
Y=taylor_documented(x,tmin,tmax,m); %running the numerical integrator

%%%%%%%%%%%%%%%%%%%%%%%%Plotting the 3D orbit
figure(1) 
plot3(Y(:,1),Y(:,2),Y(:,3))
