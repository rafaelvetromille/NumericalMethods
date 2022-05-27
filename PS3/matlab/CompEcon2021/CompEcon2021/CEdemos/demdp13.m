%% DEMDP13 Inventory Management Model
%
% Profit maximizing entrepeneur must decide how much to produce and how 
% much inventory to hold.

% States
%     s1      market price
%     s2      beginning inventory
% Actions
%     x1      quantity produced
%     x2      ending inventory
% Parameters
%     c       production cost function parameters
%     k       inventory holding cost function parameters
%     pbar	  long-run mean price
%     rho     mean-reversion coefficient
%     sigma   standard deviation of price shocks
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)



%% FORMULATION

% Numerical control parameters
basistype = 'spli';                     % basis function type
n    = [15 60];                         % degree of approximation
smin = [0.2 0.0];                     	% minimum states
smax = [2.5 3.0];                     	% maximum states
m    =     3;                           % number of shocks
nper =    26;                           % number of periods simulated
nrep = 5000;                            % number of replications
nplot =   15;                           % number of parameter values plotted in sensitivity analysis
 
% Base case model parameters
c     = [0.5 0.1];                      % production cost function parameters
k     = [0.1 0.1];                      % inventory holding cost function parameters
pbar  = 1.0;                            % long-run mean price
rho   = 0.5;                            % mean-reversion coefficient
sigma = 0.2;                            % standard deviation of price shocks
delta = 0.9;                            % discount factor

% Approximation structure
basis = fundefn(basistype,n,smin,smax);

% Continuous state shock distribution
[e,w] = qnwnorm(m,0,sigma^2);           % shocks and probabilities

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {c k pbar rho};          % function parameters
model.discount = delta;                 % discount factor
model.ds = 2;                           % dimension of continuous state
model.dx = 2;                           % dimension of continuous action
model.e  = e;                           % continuous state shocks
model.w  = w;                           % continuous state shock probabilities


%% BASE CASE MODEL SOLUTION

% Nonstochastic steady-state
sstar = [pbar 0];                       % state
xstar = [(pbar-c(1))/c(2) 0];           % action

% Check model derivatives
dpcheck(model,sstar,xstar)

% Solve collocation equation
[v,x,c,sr,vr,xr,resid] = dpsolve(model,basis);

% Reshape output for plotting
n = basis.n*10+1;
s1 = reshape(sr(:,1),n);
s2 = reshape(sr(:,2),n);
vr = reshape(vr,n);
x1 = reshape(xr(:,1),n);
x2 = reshape(xr(:,2),n);
resid = reshape(resid,n);

% Shadow prices
p1 = funeval(c,basis,sr,[1 0]);
p1 = reshape(p1,n);

% Plot optimal policy 1
figure
mesh(s1,s2,x1)
title('Optimal Production Policy')
xlabel('Market Price')
ylabel('Beginning Inventory')
zlabel('Production')

% Plot optimal policy 2
figure
mesh(s1,s2,x2)
title('Optimal Inventory Policy')
xlabel('Market Price')
ylabel('Beginning Inventory')
zlabel('Ending Inventory')

% Plot value function
figure
mesh(s1,s2,vr)
title('Value Function')
xlabel('Market Price')
ylabel('Beginning Inventory')
zlabel('Value of the Firm')

% Plot shadow price function 1
figure
mesh(s1,s2,p1)
title('Shadow Price of Inventories')
xlabel('Market Price')
ylabel('Beginning Inventory')
zlabel('Price')

% Plot residual
figure
mesh(s1,s2,resid)
title('Bellman Equation Residual')
xlabel('Market Price')
ylabel('Beginning Inventory')
zlabel('Residual')


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randnorm(0,sigma^2,nrep,nper);

% Initialize simulation
sinit = [pbar*ones(nrep,1) zeros(nrep,1)];

% Simulate model
[ssim,xsim] = dpsimul(model,basis,nper,sinit,[],sr,vr,xr,esim);
s1sim = ssim(:,:,1);
s2sim = ssim(:,:,2);
x1sim = xsim(:,:,1);

% Ergodic moments
s1avg = mean(s1sim(:)); 
s2avg = mean(s2sim(:)); 
x1avg = mean(x1sim(:)); 
s1std = std(s1sim(:)); 
s2std = std(s2sim(:)); 
x1std = std(x1sim(:)); 

% Print ergodic moments
fprintf('\n')
fprintf('Ergodic Moments\n') 
fprintf('              Nonstochastic    Ergodic      Ergodic\n') 
fprintf('              Steady-State      Mean     Std Deviation\n') 
fprintf('Market Price     %5.3f         %5.3f         %5.3f\n',[sstar(1) s1avg s1std])
fprintf('Inventory        %5.3f         %5.3f         %5.3f\n',[sstar(2) s2avg s2std])
fprintf('Production       %5.3f         %5.3f         %5.3f\n\n',[xstar(1) x1avg x1std])

% Plot simulated and expected state paths 1
figure
hold on
plot(0:nper-1,s1sim(1:3,:))
plot(0:nper-1,mean(s1sim),'k')
ytickformat('%.1f')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Market Price')

% Plot simulated and expected state paths 2
figure
hold on
plot(0:nper-1,s2sim(1:3,:))
plot(0:nper-1,mean(s2sim),'k')
ytickformat('%.1f')
title('Simulated and Expected Ending Inventory')
xlabel('Period')
ylabel('Ending Inventory')

% Plot simulated and expected action paths
figure
hold on
plot(0:nper-1,squeeze(x1sim(1:3,:)))
plot(0:nper-1,mean(x1sim),'k')
title('Simulated and Expected Production')
xlabel('Period')
ylabel('Production')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
optset('dpsolve','nr',0); 
optset('dpsolve','output',0);
s2avgplot = zeros(nplot,1);
s2stdplot = zeros(nplot,1);
vinit = v;
xinit = x;

fprintf('Varying Standard Deviation of Price Shocks\n')
modelplot = model;
sigmamin = 0.01;
sigmamax = 0.4;
sigmaplot = nodeunif(nplot,sigmamin,sigmamax);
v = vinit; x = xinit;
for ip=1:nplot
  [e,w] = qnwnorm(m,0,sigmaplot(ip)^2);
  modelplot.e = e;
  modelplot.w = w;
  rng('default')
  esim = randnorm(0,sigmaplot(ip)^2,nrep,nper);
  [v,x,c,sr,vr,xr] = dpsolve(modelplot,basis,v,x);
  [ssim,xsim] = dpsimul(modelplot,basis,nper,sinit,[],sr,vr,xr,esim);
  s2sim = ssim(:,:,2);
  s2avgplot(ip) = mean(s2sim(:));
  s2stdplot(ip) = std(s2sim(:));
end
figure
plot(sigmaplot,s2avgplot)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Standard Deviation of Price Shocks')
ylabel('Ergodic Mean Inventory')
figure
plot(sigmaplot,s2stdplot)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Standard Deviation of Price Shocks')
ylabel('Ergodic Standard Deviation of Inventory')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,x,~,~,e,c,k,pbar,rho)
n = size(s,1);
ds = 2;
dx = 2;
switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = inf*ones(size(s));
    out3 = [];
  case 'f'      % reward
    out2 = zeros(n,dx);
    out3 = zeros(n,dx,dx);
    out1 = s(:,1).*(s(:,2)+x(:,1)-x(:,2)) ...
      - (c(1)+0.5*c(2)*x(:,1)).*x(:,1) ...
      - (k(1)+0.5*k(2)*x(:,2)).*x(:,2);
    out2(:,1) =  s(:,1) - (c(1)+c(2)*x(:,1));
    out2(:,2) = -s(:,1) - (k(1)+k(2)*x(:,2));
    out3(:,1,1) = -c(2)*ones(n,1);
    out3(:,2,2) = -k(2)*ones(n,1);
  case 'g'      % transition
    out2 = zeros(n,ds,dx);
    out3 = zeros(n,ds,dx,dx);
    out1 = [pbar+rho*(s(:,1)-pbar)+e x(:,2)];
    out2(:,2,2) = ones(n,1);
end
end