%Function that implements a integrator based on a Taylor method with
%variable step h and order p
%Author: Jaime Burgos-Garcia. FCFM. Autonomous University of Coahuila
function w=taylor_documented(x,tmin,tmax,m)
t=tmin;
tolrel=2.22045e-014;
tolabs=eps;  %tolerances
lambda1=(3/2)*(1-sqrt(1-3*m+3*m^2)); %evaluation of the parameters lambda
lambda2=(3/2)*(1+sqrt(1-3*m+3*m^2));
numvar=7; %number of variables in the augmented vector field
w(1, :)=x;%saving the vector of initial conditions in a matrix w
l=2;  %index for the matrix w of approximations of the solution
i=1;  %index for the vector Times to save the steps in time 
flag=1;  %conditional

while flag==1 %main loop to compute the approximation in each step of time
%estimating the order in every step
if (tolrel*norm(w(l-1,:),Inf)<=tolabs)
   tol=tolabs;
else
    tol=tolrel;
end
p=ceil(1-(1/2)*log(tol));
    
 %%%%%% vectors to store the approximation of order p in each step of time
x1=zeros(1,p);
x2=zeros(1,p);
x3=zeros(1,p);
x4=zeros(1,p);
x5=zeros(1,p);
x6=zeros(1,p);
x7=zeros(1,p);
%%%%%% storing the first elements of the approximation
x1(1)=w(i,1);
x2(1)=w(i,2);
x3(1)=w(i,3);
x4(1)=w(i,4);
x5(1)=w(i,5);
x6(1)=w(i,6);
x7(1)=w(i,7);
 %%%%%% vectors to store the values of the auxilar variables (products) for
 %%%%%% the equations of motion
u1=zeros(1,p);
u2=zeros(1,p);
u3=zeros(1,p);
u4=zeros(1,p);
u5=zeros(1,p);
u6=zeros(1,p);
u7=zeros(1,p);
u8=zeros(1,p);
u9=zeros(1,p);
u10=zeros(1,p);
u24=zeros(1,p);
u25=zeros(1,p);
u26=zeros(1,p);
u27=zeros(1,p);
u28=zeros(1,p);
u29=zeros(1,p);
u30=zeros(1,p);
u31=zeros(1,p);

 for n=1:p %loop to compute the value of the coefficients of the partial sum of order p
   %loops to compute the values of the auxilar variables (products) for
 %%%%%% the equations of motion
    sum=0;
    for r=1:n
     sum=sum+x1(n+1-r)*x7(r);
    end
    u1(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x2(n+1-r)*x7(r);
     end
    u2(n)=sum;
    
     sum=0;
    for r=1:n
     sum=sum+x3(n+1-r)*x7(r);
     end
    u3(n)=sum;
    
     sum=0;
    for r=1:n
     sum=sum+x7(n+1-r)*x7(r);
     end
    u4(n)=sum;
 
   
      sum=0;
    for r=1:n
     sum=sum+u1(n+1-r)*u1(r);
     end
    u5(n)=sum;
    
      sum=0;
    for r=1:n
     sum=sum+u2(n+1-r)*u2(r);
     end
    u6(n)=sum;
    
     
    sum=0;
    for r=1:n
     sum=sum+u3(n+1-r)*u3(r);
     end
    u7(n)=sum;
  
    sum=0;
    for r=1:n
     sum=sum+u1(n+1-r)*u4(r);
     end
    u8(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+u2(n+1-r)*u4(r);
     end
    u9(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+u3(n+1-r)*u4(r);
     end
    u10(n)=sum;
  
    sum=0;
    for r=1:n
     sum=sum+x1(n+1-r)*x4(r);
     end
    u24(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x2(n+1-r)*x5(r);
     end
    u25(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x3(n+1-r)*x6(r);
     end
    u26(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x7(n+1-r)*u24(r);
     end
    u27(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x7(n+1-r)*u25(r);
     end
    u28(n)=sum;
    
    sum=0;
    for r=1:n
     sum=sum+x7(n+1-r)*u26(r);
     end
    u29(n)=sum;
    
   u30(n)=u27(n)+u28(n)+u29(n);
    
    sum=0;
    for r=1:n
     sum=sum+u4(n+1-r)*u30(r);
     end
    u31(n)=sum;
 %%%%%%computing the value of the coefficients 
   
 x1(n+1)=x4(n)/(n);
 x2(n+1)=x5(n)/(n);
 x3(n+1)=x6(n)/(n);
 x4(n+1)=(2*x5(n)+lambda2*x1(n)-u8(n))/(n);
 x5(n+1)=(-2*x4(n)+lambda1*x2(n)-u9(n))/(n);
 x6(n+1)=(-x3(n)-u10(n))/(n);
 x7(n+1)=-u31(n)/(n);

 end    %end of computation of the coefficients
    
   U=[x1;x2;x3;x4;x5;x6;x7]; %storing the values of the coefficients
    %estimation of the step size h
    if (tolrel*norm(w(l-1,:),Inf)<tolabs)
   rhopmenos1=(1/norm(U(:,p-1),Inf))^(1/(p-1));
else
    rhopmenos1=(norm(w(l-1,:),Inf)/norm(U(:,p-1),Inf))^(1/(p-1));
    end
    
     if (tolrel*norm(w(l-1,:),Inf)<tolabs)
   rhop=(1/norm(U(:,p),Inf))^(1/p);
else
    rhop=(norm(w(l-1,:),Inf)/norm(U(:,p),Inf))^(1/p);
     end
    rho=min([rhopmenos1 rhop]);
  
    h=(rho/exp(1)^2)*exp(-0.5/(p-1)); %choosing the step size h 
    h=min([h tmax-t]);
  
   wi=zeros(1,numvar)';
   
    for j=1:p+1 %loop to compute the approximation in current step of time
       wi=wi+U(:,j)*h^(j-1);
    end
    
     w(l, :)=wi'; %storing the computed approximation
    
        t=t+h;  %going forward to the next step of time
  l=l+1;
  i=i+1;
  
  if t>=tmax %conditional to stop the main loop
   flag=0;
  end
  
end

end
