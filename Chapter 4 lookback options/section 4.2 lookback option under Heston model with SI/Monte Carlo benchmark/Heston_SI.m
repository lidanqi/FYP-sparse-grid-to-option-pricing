function [S]=Heston_SI(dt)
Otype=2;
if Otype==1
    S=80;
else 
    S=120;
end

% define parameters
T=1;
r=0.04;
v=0.04;
rho = -0.5;

theta_v=0.02;
kappa_v=1.5;
sigma_v=0.15;
lambda_v=0;

theta_r=0.04;
kappa_r=0.3;
sigma_r=0.1;
lambda_r=0;

kappa_s=kappa_v+lambda_v;
theta_s=kappa_v*theta_v/(kappa_v+lambda_v);
kappa_r2=kappa_r+lambda_r;
theta_r2=kappa_r*theta_r/(kappa_r+lambda_r);
kappa_r=kappa_r2;
theta_r=theta_r2;


N  = floor(T/dt);
dt = T/N;

    e1 = randn(N+1, 1);
    e2 = rho*e1 + sqrt(1 - rho^2)*randn(N+1, 1); 
    e3 = randn(N+1, 1);
    
S_m=zeros(N+1,1);
v_m=zeros(N+1,1);
r_m=zeros(N+1,1);
S_m(1,:)=S;
v_m(1,:)=v;
r_m(1,:)=r; 

for i = 1:N
   
%     lnS(t) = lnS(t-1) + (mu - 0.5*V(t-1))*dt + sqrt(dt*V(t-1))*e1(t);    
%     V(t)   = max(1e-20, kappa*theta*dt + (1 - kappa*dt)*V(t-1) + sigma*sqrt(V(t-1)*dt)*e2(t));
     S_m(i+1,:)=S_m(i,:).*exp((r_m(i,:)-0.5*max(v_m(i,:),0))*dt+sqrt(max(v_m(i,:),0)).*sqrt(dt).*e1(i));
    v_m(i+1,:)=v_m(i,:)+kappa_s*(theta_s-max(v_m(i,:),0))*dt+sigma_v*sqrt(max(v_m(i,:),0)).*sqrt(dt).*e2(i);
    r_m(i+1,:)=r_m(i,:)+kappa_r*(theta_r-max(r_m(i,:),0))*dt+sigma_r*sqrt(max(r_m(i,:),0)).*sqrt(dt).*e3(i);
end 

S = S_m;