%% DEMSOC05 Stochastic Optimal Fish Harvest Model
%
% Profit maximizing fisheries owner must decide harvesting effort to exert.

% State
%     s       stock of fish
% Control
%     h       harvest effort
% parameters
%     alpha   biological growth function scale factor
%     sigma   biological growth volatility
%     H       maximum harvest effort
%     p       market price of fish
%     k       cost function parameter
%     rho     continuous discount rate

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';     % basis function type
n = 1000;               % degree of approximation
smin = 0.01;            % minimum state
smax = 1.5;             % maximum state
sinit = smax;           % initial stock of fish
T     = 25;             % time horizon
nrep  = 1000;           % number of replications

% Base case model parameters
alpha = 0.5;            % biological growth function scale factor
sigma = 0.5;            % biological growth volatility
H     = 1;              % maximum harvest effort
p     = 1;              % market price of fish
k     = 0.25;           % cost function parameter
rho   = 0.1;            % continuous discount rate

% Model structure
clear model
model.func   = @func;                           % model function file
model.rho    = rho;                             % continuous discount rate
model.params = {alpha sigma H p k};             % function file parameters

% Approximation structure
basis = fundefn(basistype,n,smin,smax);         % basis functions


%% BASE CASE MODEL SOLUTION

% Solve HJB equation by collocation
[c,s,v,x,resid] = socsolve(model,basis);

% Optimal switch point
Vs = funeval(c,basis,s,1);
sstar = interp1(Vs-p+k./s,s,0);

% plot optimal policy
figure
hold on
plot(s,x)
plottext(sstar+0.01,[],'$s^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Fish Harvest policy')
xlabel('Fish Stock')
ylabel('Rate of Harvest')

% plot value function
figure
hold on
plot(s,v)
xtickformat('%.1f')
ytickformat('%.1f')
title('Value Function')
xlabel('Fish Stock')
ylabel('Value of Operation')

% plot shadow price function
figure
hold on
plot(s,Vs)
xtickformat('%.1f')
title('Shadow Price Function')
xlabel('Fish Stock')
ylabel('Shadow price')

% plot residual
figure
hold on
plot(s,resid)
plothdash([],0)
plotvdash(sstar,0)
plottext(sstar+0.01,[],'$s^*$')
xtickformat('%.1f')
title('HJB Equation Residual')
xlabel('Fish Stock')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Simulate model
[t,ssim,~,smean,xmean] = socsimul(c,model,basis,sinit,T,nrep);

% plot simulated and expected state path
figure
plot(t,ssim,t,smean,'k')
ytickformat('%.1f')
title('Simulated and Expected Fish Stock')
xlabel('Time')
ylabel('Fish Stock')

% plot expected control path
figure
plot(t,xmean,'k')
ytickformat('%.1f')
title('Expected Rate of Harvest')
xlabel('Time')
ylabel('Rate of Harvest')

% Natural ergodic distribution
cv0 = funbase(basis,s)\((20*p)*s);
[cp0,Ex0] = itodensity(model,basis,cv0);
p0 = funeval(cp0,basis,s);

% Ergodic distribution with harvesting
[cp,Ex] = itodensity(model,basis,c);
p = funeval(cp,basis,s);

% plot ergodic distributions
figure
hold on
plot(s,p,s,p0)
plotvdash(sstar,funeval(cp,basis,sstar))
plottext(sstar+0.01,[],'$s^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Ergodic Distribution of Fish Stock')
xlabel('Fish Stock')
ylabel('Probability')
legend('Harvested','Natural')

% print switching point and ergodic moments
fprintf('\n')
fprintf('Switching point and Ergodic Moments\n')
fprintf('  s*                 %5.3f\n',sstar)
fprintf('  E[s]               %5.3f\n',Ex)
fprintf('  E[s|x=0]           %5.3f\n',Ex0)
fprintf('  pct time inactive  %5.3f\n',funeval(cp,basis,sstar,-1))


%% SAVE FIGURES
printfigures(mfilename)


%% SOCSOLVE FUNCTION FILE
function out = func(flag,s,x,Vs,~,alpha,sigma,H,p,k)
switch flag
  case 'x'        % optimal control    
    out = H*(Vs<(p-k./s));
  case 'f'        % reward
    out = (p-k./s).*s.*x;
  case 'mu'       % state drift
    out = (alpha*(1-s)-x).*s;
  case 'sigma'    % state diffusion
    out = sigma*s;
end
end