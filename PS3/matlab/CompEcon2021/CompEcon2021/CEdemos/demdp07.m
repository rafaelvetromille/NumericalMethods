%% DEMDP07 Solow Stochastic Optimal Economic Growth Model
%
% Welfare maximizing social planner must decide how much society should
% consume and invest. Unlike the Ramsey model, this model allows arbitrary
% constant relative risk aversion and capital depreciation.  Unlike the
% Ramsey model, it lacks a known closed-form solution.

% States
%     s       stock of wealth
% Actions
%     k       capital investment
% Parameters
%     alpha   relative risk aversion
%     beta    production elasticity
%     gamma   capital survival rate
%     sigma   production shock volatility
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';                     % basis function type
n    =    21;                           % degree of approximation
smin =     5;                           % minimum stock of wealth
smax =    20;                           % maximum stock of wealth
m    =     5;                           % number of production shocks
nper =    81;                           % number of periods simulated
nrep =  5000;                           % number of replications
nplot =   25;                           % number of parameter values plotted
    
% Base case model parameters
alpha = 2.0;                            % relative risk aversion
beta  = 0.5;                            % production elasticity
gamma = 0.9;                          	% capital survival rate
sigma = 0.1;                            % production shock volatility
delta = 0.9;                            % discount factor

% Approximation structure
[basis,~,snodes] = fundefn(basistype,n,smin,smax);

% Production shock distribution
[e,w] = qnwlogn(m,1,sigma^2);

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {alpha beta gamma};	    % function parameters
model.discount = delta;                 % discount factor
model.e  = e;                           % shocks
model.w  = w;                           % shock probabilities


%% BASE CASE MODEL SOLUTION
  
% Nonstochastic steady-state
kstar = ((1-delta*gamma)/(delta*beta))^(1/(beta-1));  % capital investment
sstar = gamma*kstar + kstar^beta;                     % wealth

% Check model derivatives
dpcheck(model,sstar,kstar)

% Initialize v and k via LQ approximation
[vlq,klq] = lqapprox(model,snodes,smin,smin*kstar/sstar);

% Solve collocation equation
[v,k,c,sr,vr,kr,resid] = dpsolve(model,basis,vlq,klq);

% Plot optimal policy
figure
hold on
plot(sr,kr)
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Investment')

% Plot value function
figure
plot(sr,vr)
ytickformat('%.1f')
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')

% Plot shadow price function
figure
pr = funeval(c,basis,sr,1);
plot(sr,pr)
ytickformat('%.1f')
title('Shadow Price Function')
xlabel('Wealth')
ylabel('Shadow Price')

% Plot collocation equation residual
figure
hold on
plot(sr,resid)
plothdash([],0)
xlim([smin smax])
ytickformat('%.1f')
title('Collocation Equation Residual')
xlabel('Wealth')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(1,sigma^2,nrep,nper);

% Initialize simulation
sinit = smin*ones(nrep,1);

% Simulate model
[ssim,ksim] = dpsimul(model,basis,nper,sinit,[],sr,vr,kr,esim);

% Plot simulated and expected wealth path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
ytickformat('%.1f')
title('Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

% Plot simulated and expected investment path
figure
hold on
plot(0:nper-1,ksim(1:3,:))
plot(0:nper-1,mean(ksim),'k')
ytickformat('%.1f')
title('Simulated and Expected Investment')
xlabel('Period')
ylabel('Investment')

% Ergodic moments
ssim = ssim(:,nper);
ksim = ksim(:,nper);
savg = mean(ssim); 
kavg = mean(ksim); 
sstd = std(ssim); 
kstd = std(ksim); 

% Print stead-state and ergodic moments
fprintf('          Nonstochastic    Ergodic      Ergodic\n') 
fprintf('          Steady-State      Mean     Std Deviation\n') 
fprintf('Wealth       %5.3f         %5.3f         %5.3f\n'  ,[sstar savg sstd])
fprintf('Investment   %5.3f         %5.3f         %5.3f\n\n',[kstar kavg kstd])

% Plot ergodic wealth distribution
figure
histogram(ssim(:),'Normalization','pdf')
xlim([5 10])  
ytickformat('%.1f')
title('Ergodic Wealth Distribution')
xlabel('Wealth')
ylabel('Probability Density')

% Plot kernel-smoothed ergodic wealth distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
xlim([5 10])  
ytickformat('%.1f')
title('Ergodic Wealth Distribution')
xlabel('Wealth')
ylabel('Probability Density')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initializations
savgplot = zeros(nplot,1);
sstdplot = zeros(nplot,1);
vinit = v;
kinit = k;

fprintf('Varying Relative Risk Aversion\n')
modelplot = model;
alphamin = 0.3;
alphamax = 4.0;
alphaplot = nodeunif(nplot,alphamin,alphamax);
v = vinit; k = kinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {alphaplot(ip) beta gamma};
  [v,k,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(alphaplot,savgplot)
ytickformat('%.2f')
xlabel('Relative Risk Aversion')
ylabel('Ergodic Mean Wealth')
figure
plot(alphaplot,sstdplot)
ytickformat('%.1f')
xlabel('Relative Risk Aversion')
ylabel('Ergodic Standard Deviation of Wealth')

fprintf('Varying Production Elasticity\n')
modelplot = model;
betamin = 0.3;
betamax = 0.6;
betaplot = nodeunif(nplot,betamin,betamax);
v = vinit; k = kinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {alpha betaplot(ip) gamma};
  [v,k,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(betaplot,savgplot)
xlabel('Production Elasticity')
ylabel('Ergodic Mean Wealth')
xticks(betamin:0.1:betamax)
figure
plot(betaplot,sstdplot)
xlabel('Capital Production Elasticity')
ylabel('Ergodic Standard Deviation of Wealth')
xticks(betamin:0.1:betamax)
ytickformat('%.1f')

fprintf('Varying Production Shock Volatility\n')
modelplot = model;
sigmamin = 0.0;
sigmamax = 0.3;
sigmaplot = nodeunif(nplot,sigmamin,sigmamax);
v = vinit; k = kinit;
for ip=1:nplot
  [e,w] = qnwlogn(m,1,sigmaplot(ip)^2);
  modelplot.e  = e;
  modelplot.w  = w;
  rng('default')
  esim = randlogn(1,sigmaplot(ip)^2,nrep,nper);
  [v,k,c,sr,vr,kr] = dpsolve(modelplot,basis,v,k);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,kr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(sigmaplot,savgplot)
xlabel('Production Shock Volatility')
ylabel('Ergodic Mean Wealth')
xticks(sigmamin:0.1:sigmamax)
ytickformat('%.1f')
figure
plot(sigmaplot,sstdplot)
xlabel('Production Shock Volatility')
ylabel('Ergodic Standard Deviation of Wealth')
xticks(sigmamin:0.1:sigmamax)
ytickformat('%.1f')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,k,~,~,e,alpha,beta,gamma)
switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = 0.99*s;
    out3 = [];
  case 'f'      % reward
    out1 = ((s-k).^(1-alpha))/(1-alpha);
    out2 = -(s-k).^(-alpha);
    out3 = -alpha*(s-k).^(-alpha-1);
  case 'g'      % transition
    out1 = gamma*k + e.*k.^beta;
    out2 = gamma + beta*e.*k.^(beta-1);
    out3 = (beta-1)*beta*e.*k.^(beta-2);
end
end