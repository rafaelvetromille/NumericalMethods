%% DEMRS01 Asset Abandonment Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION
  
% Base case model parameters
c     = 0.5;
mu    = 0;
sigma = 0.2;
rho   = 0.1;

% Model structure
clear model
model.func   = @func;
model.params = {c mu sigma rho};
model.xindex = [1 0 1 2 0;1 0 0 0 1];


%% BASE CASE MODEL SOLUTION

% Initial values and approximation size
n = 101;
x = [0;20];

% Solve collocation equation
[cv,basis,x] = rssolve(model,x,n,'cheb');
cv = cv{1};
basis = basis{1};
sstar = x(1);

% Plot results
S = nodeunif(1001,0,1);
V   = funeval(cv,basis,S);
Vs  = funeval(cv,basis,S,1);
Vss = funeval(cv,basis,S,2);
V(S<sstar)   = 0;
Vs(S<sstar)  = 0;
Vss(S<sstar) = 0;
vstar  = funeval(cv,basis,sstar);
vstar1 = funeval(cv,basis,sstar,1);

% Plot value function
figure
hold on
plot(S,V)
plotvdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$p^*$')
xtickformat('%.1f')
title('Value Function')
xlabel('$P$')
ylabel('$V$')
ylim([0 6])

% Plot first derivative of value function
figure
hold on
plot(S,Vs)
plotvdash(sstar,vstar1)
plotbullet(sstar,vstar1)
plottext(sstar+0.01,[],'$p^*$')
xtickformat('%.1f')
title('First Derivative of Value Function')
xlabel('$P$')
ylabel('$V''$')

% Exact solution
beta = roots([sigma^2/2 mu-sigma^2/2 -rho]);
beta = beta(beta<0);
pstar = (rho-mu)*beta*c/rho/(beta-1);
A = -pstar.^(1-beta)/(rho-mu)/beta;
VV = S/(rho-mu)-c/rho+A*S.^beta;
VV(S<pstar) = 0;

% Exact approximation error
e = rho*V-(S-c).*(S>sstar)-mu*S.*Vs-0.5*sigma^2*S.^2.*Vss;

% Plot residual and approximation error
figure
hold on
plot(S,e)
plot(S,VV-V(:,end))
plothdash([],0)
xtickformat('%.1f')
ytickformat('%.2f')
title('Approximation Residual')
xlabel('$P$')
ylabel('Residual/Residual')
legend('Residual','Approximation Error')

fprintf('\n')
fprintf('Switching Point\n') 
fprintf('      "Exact"      Approximate        Error\n') 
fprintf('       %5.3f         %5.3f         %8.1e\n',[sstar(1) pstar(1) sstar(1)-pstar(1)])


%% SAVE FIGURES
printfigures(mfilename)


%% RSSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,x,<params>)

function out = func(flag,s,~,c,mu,sigma,rho)
switch flag
  case 'f'        % reward
    out = (s-c);
  case 'g'        % state drift
    out = mu*s;
  case 'sigma'    % state diffusion
    out = sigma*s;
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'reward'   % jump rewards and derivatives
    out = [0 0 0;s(2)/(rho-mu)-c/rho 1/(rho-mu) 0];
end
end