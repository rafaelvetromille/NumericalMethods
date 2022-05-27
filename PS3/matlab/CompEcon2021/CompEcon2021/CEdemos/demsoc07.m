%% DEMSOC07 Stochastic Nonrenewable Resource Model
%
% Welfare maximizing social planner must decide the rate at which a
% nonrenewable resource should be harvested.

% State
%     s       resource stock
% Control
%     q       harvest rate
% Parameters
%     kappa   harvest unit cost scale factor
%     gamma   harvest unit cost elasticity
%     eta     inverse elasticity of demand
%     rho     discount rate
%     sigma   diffusion volatility

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';     % basis function type
n     = 31;             % degree of approximation
smin  = 0.1;            % minimum resource stock
smax  = 1.0;            % maximum resource stock
sinit = smax;           % initial capital stock
T     = 20;             % time horizon
nrep  = 1000;           % number of replications

% Base case model parameters
kappa = 10;             % harvest unit cost scale factor
gamma = 1;              % harvest unit cost elasticity
eta   = 1.5;            % inverse elasticity of demand
rho   = 0.05;           % discount rate
sigma = 0.1;            % diffusion volatiity

% Model Structure
clear model
model.func   = @func;                           % model function file
model.rho    = rho;                             % discount rate
model.params = {eta kappa gamma sigma};         % function file parameters

% Approximation Structure
basis = fundefn(basistype,n,smin,smax);         % basis functions


%% BASE CASE MODEL SOLUTION

% Solve HJB equation by collocation
[c,s,v,q,resid] = socsolve(model,basis);

% Plot optimal policy
figure
plot(s,q)
xtickformat('%.1f')
ytickformat('%.3f')
title('Optimal Harvest Policy')
xlabel('Resource Stock')
ylabel('Rate of Harvest')

% Plot value function
figure
plot(s,v)
xlim([0 smax])
xtickformat('%.1f')
title('Value Function')
xlabel('Resource Stock')
ylabel('Social Welfare')

% Plot shadow price function
figure
Vs = funeval(c,basis,s,1);
plot(s,Vs)
xtickformat('%.1f')
title('Shadow Price Function')
xlabel('Resource Stock')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
xtickformat('%.1f')
title('HJB Equation Residual')
xlabel('Resource Stock')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Simulate model
[t,ssim,xsim,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);

% Plot simulated and expected state path
figure
plot(t,ssim,t,smean,'k')
ytickformat('%.1f')
title('Simulated and Expected Resource Stock')
xlabel('Time')
ylabel('Resource Stock')

% Plot simulated and expected control path
figure
plot(t,xsim,t,xmean,'k')
ytickformat('%.3f')
title('Simulated and Expected Rate of Harvest')
xlabel('Time')
ylabel('Rate of Harvest')


%% SAVE FIGURES
printfigures(mfilename)


%% SOCSOLVE FUNCTION FILE
function out = func(flag,s,q,Vs,~,eta,kappa,gamma,sigma)
switch flag
  case 'x'        % optimal control
    k = kappa*s.^(-gamma);
    out = (Vs+k).^(-1/eta);
  case 'f'        % reward
    u = (1/(1-eta))*q.^(1-eta);
    k = kappa*s.^(-gamma);
    out = u - k.*q;
  case 'mu'       % state drift
    out = -q;
  case 'sigma'    % state diffusion
    out = sigma*s;
end
end