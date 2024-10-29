function du = lorenzfield_normalized(t,u)
s=10;
r=28;
b=8/3;
norm=sqrt((s*(u(2)-u(1)))^2+(r*u(1)-u(1)*u(3)-u(2))^2+(u(1)*u(2)-b*u(3))^2);
du=[(s*(u(2)-u(1)))/norm; (r*u(1)-u(1)*u(3)-u(2))/norm; (u(1)*u(2)-b*u(3))/norm; 1/norm];

end