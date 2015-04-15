% MC pricing of the European Call
function result=MCpricingHestLB(Num_simulation)
tic
% fprintf(['lookback call, fCIR, many_results, p59\n' ...
%          'no variance control\n']);


    S=10;

Mm=[6,7,8,9,10];

K    = 10;
T    = 1;
rate = 0.05;

a = 0.1;
b = 0.5;
y = 0.25;
rho   = 0;

S0 = 10;
V0 = 0.5;
params = [a, b, y, rho];


% call1 = Heston_Call(params, S0, K, T, rate, V0);
call1=4.6;

timer = clock;
payoff = zeros(Num_simulation, 5);
payoff_var = zeros(Num_simulation,5);
% rng(314159);
% 
%     [S] = Heston_Sim_Euler([rate,params], S0, V0, T, 2^(-13),'setSeed',1);
%     payoff(1) = max(S(end)-min(S);
parfor i = 1:Num_simulation
    [S] = Heston_Sim_Euler([rate,params], S0, V0, T, 2^(-13));
    Smin=min(S);
    payoff_temp = [max(S(end)-min(Smin,6), 0),max(S(end)-min(Smin,7), 0),max(S(end)-min(Smin,8), 0),max(S(end)-min(Smin,9), 0),max(S(end)-min(Smin,10), 0)];
    payoff(i,:) = payoff_temp;
    payoff_var(i,:)=max(0,S(end)-Mm);
end
% disc_payoff = exp(-rate*T)*payoff;
% call2 = mean(disc_payoff,1);
% error = std(disc_payoff)/sqrt(Num_simulation);
% est_error = abs(call1-call2);
% timespent = etime(clock,timer);
% fprintf('Estimation  : '); fprintf('%.4f  ',call2);% fprintf('Stand Error : %.6f, est error: %.6f\n',error, est_error);
% fprintf('            ( %.6f, %.6f )\n',call2-1.96*error,call2+1.96*error);
% fprintf('\ntime spent  : %.2f s\n',timespent);

disc_payoff = exp(-rate*T)*payoff;
price = mean(disc_payoff,1);
disc_payoff_var = exp(-rate*T)*payoff_var;
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

fprintf('\n%g paths, %g timesteps. \nTime %g s\n',Num_simulation,2^13,toc);


result = [price];

