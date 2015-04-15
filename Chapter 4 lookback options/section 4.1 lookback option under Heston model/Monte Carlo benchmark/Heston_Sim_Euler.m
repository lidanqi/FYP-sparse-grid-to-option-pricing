function [S, V] = Heston_Sim_Euler(params, S0, V0, T, dt,varargin)
 p = inputParser;
 p.addParamValue('setSeed', 0);
 p.parse(varargin{:});
 setSeed=p.Results.setSeed;
 
mu    = params(1);
kappa = params(2);
theta = params(3);
sigma = params(4);
rho   = params(5);

N  = floor(T/dt);
dt = T/N;

V   = zeros(N+1, 1);
lnS = zeros(N+1, 1);

V(1)   = V0;
lnS(1) = log(S0);

    w = randn(N+1, 1);
    z = rho*w + sqrt(1 - rho^2)*randn(N+1, 1);    


for t = 2:N+1
    
    lnS(t) = lnS(t-1) + (mu - 0.5*V(t-1))*dt + sqrt(dt*V(t-1))*w(t);    
    V(t)   = max(1e-20, kappa*theta*dt + (1 - kappa*dt)*V(t-1) + sigma*sqrt(V(t-1)*dt)*z(t));
    
end 

S = exp(lnS);