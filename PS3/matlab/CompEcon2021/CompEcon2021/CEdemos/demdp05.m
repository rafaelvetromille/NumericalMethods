%% DEMDP05 American Put Option Pricing Model
%
% Compute the critical exercise price for an American put option in terms 
% of time to expiration.

% States
%     p       underlying asset price
% Actions
%     j       exercize (2) or do not exercize (1)
% Parameters
%     K       option strike price
%     N       number of periods to expiration
%     mu      mean of log price innovation
%     sigma   asset price volatility
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Continuous time discretization
% sigma = 0.2;                            % annual volatility
% T     = 0.5;                            % years to expiration
% K     = 1.0;                            % option strike price
% r     = 0.1;                            % annual interest rate
% N     = 300;                            % number of time intervals
% dt    = T/N;                            % length of time interval
% delta = exp(-r*dt);                     % per period discount factor
% mu    = dt*(r-sigma^2/2);               % mean of log price innovation

% Numerical control parameters
basistype = 'spli';                     % basis function type
n    = 500;                             % degree of approximation
pmin = -1;                              % minimum log price
pmax =  1;                              % maximum log price
m    = 15;                              % number of shocks
  
% Base case model parameters
K     = 1.0;                            % option strike price
N     = 300;                            % number of periods to expiration
mu    = 0.0001;                         % mean of log price innovation
sigma = 0.0080;                         % asset price volatility
delta = 0.9998;                         % discount factor  

% Approximation structure
[basis,Phi,p] = fundefn(basistype,n,pmin,pmax);

% Continuous state shock distribution
[e,w] = qnwnorm(m,mu,sigma^2);          % shocks and probabilities  


%% BASE CASE MODEL SOLUTION
  
% Intialize value function approximant coefficients
c = zeros(n,N+1);

% Solve collocation equation by backward recursion
for t=N:-1:1
  v = zeros(n,1);
  for k=1:m
    pnext = p + e(k);
    v = v + w(k)*max(K-exp(pnext),delta*funeval(c(:,t+1),basis,pnext));
  end
  c(:,t) = Phi\v;
end

% Critical exercise prices
pcrit = zeros(N+1,1);
f = @(p,K,delta,c,basis) K-exp(p)-delta*funeval(c,basis,p);
for t=1:N
  pcrit(t) = broyden(f,0,K,delta,c(:,t),basis);
end

% Plot critical exercise prices
figure
time = (N:-1:0);
plot(time,exp(pcrit))
ytickformat('%.2f')
title('American Put Option Optimal Exercise Boundary')
xlabel('Periods Remaining Until Expiration')
ylabel('Exercise Price')


%% SAVE FIGURES
printfigures(mfilename)