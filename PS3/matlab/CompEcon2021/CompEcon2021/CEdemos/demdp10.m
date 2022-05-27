%% DEMDP10 Water Resource Management Model
%
% Public authority must decide how much water to release from a reservoir so 
% as to maximize benefits derived from agricultural and recreational uses.

% States
%     s       reservoir level at beginning of summer
% Actions
%     x       quantity of water released for irrigation
% Parameters
%     a       producer benefit function parameters
%     b       recreational user benefit function parameters
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';     % basis function type
n    =    15;           % degree of approximation
smin =     2;           % minimum reservoir level
smax =     8;           % maximum reservoir level
m    =     3;           % number of shocks
nper =    31;           % number of periods simulated
nrep = 10000;           % number of replications
nplot =   15;           % number of parameter values plotted

% Base case model parameters
a = [1 -2];             % producer benefit function parameters
b = [2 -3];             % recreational user benefit function parameters
ymean = 1.0;            % mean rainfall
sigma = 0.2;            % rainfall volatility
delta = 0.9;            % discount factor

% Continuous state shock distribution
[e,w] = qnwlogn(m,ymean,sigma^2);

% Model structure
clear model
model.func = @func;         % model functions
model.params = {a b};       % function parameters
model.discount = delta;     % discount factor
model.e  = e;               % continuous state shocks
model.w  = w;               % continuous state shock probabilities

% Approximation structure
[basis,~,s] = fundefn(basistype,n,smin,smax); 


%% BASE CASE MODEL SOLUTION

% Nonstochastic steady-state
xstar = 1;                                  % action
sstar = 1+(a(1)*(1-delta)/b(1))^(1/b(2));   % state

% Check model derivatives
dpcheck(model,sstar,xstar)

% Compute linear-quadratic approximation at collocation nodes
[vlq,xlq] = lqapprox(model,s,sstar,xstar); 

% Solve collocation equation
[v,x,c,sr,vr,xr,resid] = dpsolve(model,basis,vlq,xlq);

% Plot optimal policy
figure
plot(sr,xr)
ytickformat('%.1f')
title('Optimal Irrigation Policy')
xlabel('Reservoir Level')
ylabel('Irrigation')

% Plot value function
figure
plot(sr,vr)
title('Value Function')
xlabel('Reservoir Level')
ylabel('Social Welfare')

% Plot shadow price function
figure
pr = funeval(c,basis,sr,1);
plot(sr,pr)
ytickformat('%.1f')
title('Shadow Price Function')
xlabel('Reservoir Level')
ylabel('Shadow Price')

% Plot residual
figure
hold on
plot(sr,resid)
plothdash([],0)
ytickformat('%.1f')
title('Bellman Equation Residual')
xlabel('Reservoir Level')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(1,sigma^2,nrep,nper);

% Initialize simulation
sinit = smin*ones(nrep,1);              % initial reservoir levels

% Simulate model
[ssim,xsim] = dpsimul(model,basis,nper,sinit,[],sr,vr,xr,esim);

% Plot simulated state path
figure
hold on
plot(0:nper-1,ssim(1:3,:))
plot(0:nper-1,mean(ssim),'k')
title('Simulated and Expected Reservoir Level')
xlabel('Year')
ylabel('Reservoir Level')
ytickformat('%.1f')

% Plot simulated action path
figure
hold on
plot(0:nper-1,xsim(1:3,:))
plot(0:nper-1,mean(xsim),'k')
title('Simulated and Expected Irrigation')
xlabel('Year')
ylabel('Irrigation')
ytickformat('%.1f')

% Ergodic moments
ssim = ssim(:,nper);
xsim = xsim(:,nper);
savg = mean(ssim(:)); 
xavg = mean(xsim(:)); 
sstd = std(ssim(:)); 
xstd = std(xsim(:)); 

% Print ergodic moments
fprintf('\n')
fprintf('Ergodic Moments\n') 
fprintf('                 Nonstochastic    Ergodic      Ergodic\n') 
fprintf('                 Steady-State      Mean     Std Deviation\n') 
fprintf('Reservoir Level     %5.2f         %5.2f         %5.2f\n'  ,[sstar savg sstd])
fprintf('Irrigation          %5.2f         %5.2f         %5.2f\n\n',[xstar xavg xstd])

% Plot ergodic water reservoir level distribution
[qq,ss] = ksdensity(ssim(:),'support','positive');
figure
plot(ss,qq)
title('Ergodic Reservoir Level Distribution')
xlabel('Reservoir Level')
ylabel('Density')
xlim([2 6])  
ytickformat('%.1f')

% Plot ergodic water reservoir level distribution
figure
histogram(ssim(:),'Normalization','pdf')
title('Ergodic Reservoir Level Distribution')
xlabel('Reservoir Level')
ylabel('Density')
xlim([2 6])  
ytickformat('%.1f')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Intialize for sensitivity analysis
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
savgplot = zeros(nplot,1);
sstdplot = zeros(nplot,1);
vinit = v;
xinit = x;

fprintf('Varying Mean Rainfall\n')
modelplot = model;
ymeanmin = 0.8;
ymeanmax = 1.2;
ymeanplot = nodeunif(nplot,ymeanmin,ymeanmax);
v = vinit; x = xinit;
for ip=1:nplot
  [e,w] = qnwlogn(m,ymeanplot(ip),sigma^2);
  modelplot.e = e;
  modelplot.w = w;
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,v,x);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,xr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(ymeanplot,savgplot)
xticks(ymeanmin:0.1:ymeanmax)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Mean Rainfall')
ylabel('Ergodic Mean Reservoir Level')
figure
plot(ymeanplot,sstdplot)
xticks(ymeanmin:0.1:ymeanmax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Mean Rainfall')
ylabel('Ergodic Standard Deviation of Reservoir Level')

fprintf('Varying Welfare Weight on Farmers\n')
modelplot = model;
fweightmin = 0.9;
fweightmax = 1.1;
fweightplot = nodeunif(nplot,fweightmin,fweightmax);
v = vinit; x = xinit;
for ip=1:nplot
  rng('default')
  modelplot.params = {fweightplot(ip)*a b};
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,v,x);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,xr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(fweightplot,savgplot)
xlim([fweightmin fweightmax])
xticks(fweightmin:0.1:fweightmax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Welfare Weight on Farmers')
ylabel('Ergodic Mean Reservoir Level')
figure
plot(fweightplot,sstdplot)
xlim([fweightmin fweightmax])
xticks(fweightmin:0.1:fweightmax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Welfare Weight on Farmers')
ylabel('Ergodic Standard Deviation of Reservoir Level')

fprintf('Varying Rainfall Volatillity\n')
modelplot = model;
sigmamin = 0.1;
sigmamax = 0.4;
sigmaplot = nodeunif(nplot,sigmamin,sigmamax);
v = vinit; x = xinit;
for ip=1:nplot
  [e,w] = qnwlogn(m,ymean,sigmaplot(ip)^2);
  modelplot.e = e;
  modelplot.w = w;
  rng('default')
  esim = randlogn(1,sigmaplot(ip)^2,nrep,nper);
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,v,x);
  ssim = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,xr,esim);
  savgplot(ip) = mean(ssim(:,nper));
  sstdplot(ip) = std(ssim(:,nper));
end
figure
plot(sigmaplot,savgplot)
xticks(sigmamin:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Rainfall Volatillity')
ylabel('Ergodic Mean Reservoir Level')
figure
plot(sigmaplot,sstdplot)
xticks(sigmamin:0.1:sigmamax)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Rainfall Volatillity')
ylabel('Ergodic Standard Deviation of Reservoir Level')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,x,~,~,e,a,b)
switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (a(1)/(1+a(2)))*x.^(1+a(2))+(b(1)/(1+b(2)))*(s-x).^(1+b(2));
    out2 = a(1)*x.^a(2)-b(1)*(s-x).^b(2);
    out3 = a(1)*a(2)*x.^(a(2)-1)+b(1)*b(2)*(s-x).^(b(2)-1);
  case 'g'      % transition
    out1 = s-x+e;
    out2 = -ones(size(s));
    out3 = zeros(size(s));
end
end