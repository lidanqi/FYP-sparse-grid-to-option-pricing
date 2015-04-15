% [a,b]=Heston_MCS(80,100,1,0.03,0.04,2,0.04,0,0.2,0.5)
function [price] = Heston_MCS_LB2(Num_simulation)
tic
% type
Otype=2;


if Otype==1
    S0=80;
else 
    S0=120;
end
Mm=[80,90,100,110,120];

  

% define parameters
T=1;
r=0.04;
% v=0.04;

% theta_v=0.02;
% kappa_v=1.5;
% sigma_v=0.15;
% lambda_v=0;

% rho = -0.5;

% theta_r=0.04;
% kappa_r=0.3;
% sigma_r=0.1;
% lambda_r=0;


% start running

% Num_simulation=10000;
% kappa_s=kappa_v+lambda_v;
% theta_s=kappa_v*theta_v/(kappa_v+lambda_v);

% kappa_r2=kappa_r+lambda_r;
% theta_r2=kappa_r*theta_r/(kappa_r+lambda_r);
% kappa_r=kappa_r2;
% theta_r=theta_r2;

% num_Mm=size(Mm,2);

N=2^13;
payoff = zeros(Num_simulation, 5);
payoff_var = zeros(Num_simulation,5);
K=Mm;
parfor i=1:Num_simulation
    [S] = Heston_SI(1/N);
    Smin=min(S);
    payoff_temp=max(0, S(end)-min(Smin,Mm));
%     payoff_temp = [max(max(Smax,10)-S0(end),0),max(max(Smax,11)-S0(end),0),max(max(Smax,12)-S0(end),0),max(max(Smax,13)-S0(end),0),max(max(Smax,14)-S0(end),0),];
    payoff(i,:) = payoff_temp;
    payoff_var(i,:)=max(0,S(end)-K);
end

disc_payoff = exp(-r*T)*payoff;
price = mean(disc_payoff,1);
disc_payoff_var = exp(-r*T)*payoff_var;
price_var = mean(disc_payoff_var,1);

error_b4var = std(disc_payoff)/sqrt(Num_simulation);
% determine the beta:
beta=[0 0 0 0 0];
price_final=price;
disc_payoff_final=disc_payoff;
for i=1:5
beta(i)=sum((disc_payoff(:,i)-price(:,i)).*(disc_payoff_var(:,i)-price_var(:,i))) / sum( (disc_payoff_var(:,i)-price_var(:,i)).^2);
disc_payoff_final(:,i) = disc_payoff(:,i)-beta(i)*(disc_payoff_var(:,i) - price_var(:,i));
end
% estimate with formula: 

price_final=mean(disc_payoff_final,1);
sd=std(disc_payoff_final);
se2=sd/sqrt(Num_simulation);
price=price_final;

fprintf(  'm     =  ');
fprintf('%.1f  ',Mm);
fprintf('\nPrice Esti =  ');
fprintf('%.4f  ',price);
fprintf('\nError b4   =  ');
fprintf('%.4f  ',error_b4var);
fprintf('\nError var  =  ');
fprintf('%.4f  ',se2);

fprintf('\n%g paths, %g timesteps. \nTime %g s\n',Num_simulation,N,toc);

if (Num_simulation>=10000)
    load handel
    sound(y,Fs)
    if (Num_simulation>=20000)
          details=zeros(4,5);
          a=clock;
          date=a(2)*1000000+a(3)*10000+a(4)*100+a(5);
          details(1,1)=date;details(1,2)=Num_simulation;
          details(2,:)=price;
          details(3,:)=se2;
          details(4,:)=error_b4var;
          save('MCresults2','details');
    end
end
if (Num_simulation<10000 && (Num_simulation>=1000))
    load chirp
    sound(y,Fs)
end
% 
% parfor idx=1:Num_simulation
%     
% S_m=zeros(N+1,1);
% v_m=zeros(N+1,1);
% r_m=zeros(N+1,1);
% S_m(1,:)=S0;
% v_m(1,:)=v;
% r_m(1,:)=r; 
%      mat_m=nan(num_Mm,1);
%     mat_ST=nan(num_Mm,1);   
% for i=1:N    
% %     e1=norminv(random('unif',0,1),0,1);
% %     e2_temp=norminv(random('unif',0,1),0,1);
%     e1=randn(1,1);
%     e2_temp=randn(1,1);
%     e2=e1*rho+e2_temp*sqrt(1-rho*rho);
%     e3=randn(1,1);
%     
%     S_m(i+1,:)=S_m(i,:).*exp((r_m(i,:)-0.5*max(v_m(i,:),0))*dt+sqrt(max(v_m(i,:),0)).*sqrt(dt).*e1);
%     v_m(i+1,:)=v_m(i,:)+kappa_s*(theta_s-max(v_m(i,:),0))*dt+sigma_v*sqrt(max(v_m(i,:),0)).*sqrt(dt).*e2;
%     r_m(i+1,:)=r_m(i,:)+kappa_r*(theta_r-max(r_m(i,:),0))*dt+sigma_r*sqrt(max(r_m(i,:),0)).*sqrt(dt).*e3;
% end
% 
% if (Otype==1)
%     Smax=max(S_m,[],1);
%     Smax=max(Smax,Mm);
%     C=exp(-r*T)*(Smax-S_m(N+1,:));
% else
%     Smin=min(S_m,[],1);
%     for index=1:num_Mm
%      mat_m(index,:)=min(Smin,Mm(index));
%      mat_ST(index,:)=S_m(N+1,:);
%     end
% %     Smin=min(Smin,Mm);
%     C(:,idx)=exp(-r*T)*(mat_ST-mat_m);
% end 
%     C(:,idx)=exp(-r*T)*(mat_ST-mat_m);
% end
% price=mean(C,2);
% err=std(C)/sqrt(Num_simulation);
% fprintf('S0=%g, price=%.4f,\nerror: %.5f, CI=[%.4f, %.4f]\n',S,price,abs(price-bc),price-1.96*err,price+1.96*err);

% fprintf(  'm     =  ');
% fprintf('%.1f  ',Mm);
% fprintf('\nPrice =  ');
% fprintf('%.4f  ',price);
% fprintf('\n%g paths, %g timesteps. \nTime %g s\n',Num_simulation,N,toc);