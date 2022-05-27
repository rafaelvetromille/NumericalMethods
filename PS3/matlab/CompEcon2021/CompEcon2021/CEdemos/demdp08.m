%% DEMDP08 Public Renewable Resource Model
%
% Welfare maximizing social planner must decide how much of a renewable 
% resource to harvest.

% States
%     s       quantity of stock available
% Actions
%     q       quantity of stock harvested
% Parameters
%     alpha   growth function parameter
%     beta    growth function parameter
%     gamma   relative risk aversion
%     kappa   unit cost of harvest
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';   % basis function type
n     =  8;           % degree of approximation
smin  =  6;           % minimum stock
smax  =  9;           % maximum stock
nper  = 11;           % number of periods simulated
nplot = 35;           % number of parameter values plotted

% Base case model parameters
alpha = 4.0;          % growth function parameter
beta  = 1.0;          % growth function parameter
gamma = 0.5;          % relative risk aversion
kappa = 0.2;          % unit cost of harvest
delta = 0.9;          % discount factor

% Approximation structure
basis = fundefn(basistype,n,smin,smax);
    
% Model structure
clear model
model.func = @func;                             % model functions
model.params = {alpha beta gamma kappa};        % function parameters
model.discount = delta;                         % discount factor


%% BASE CASE MODEL SOLUTION
  
% Steady-state
sstar = (alpha^2-1/delta^2)/(2*beta);           % stock
qstar = sstar - (delta*alpha-1)/(delta*beta); 	% action
pstar = qstar.^(-gamma);                        % market price
lstar = pstar-kappa;                            % shadow price
vstar = ((qstar.^(1-gamma))/(1-gamma)-kappa*qstar)/(1-delta); % value
fprintf('Steady States\n') 
fprintf('   Stock         %5.2f\n'  ,sstar)
fprintf('   Harvest       %5.2f\n'  ,qstar)
fprintf('   Market Price  %5.2f\n'  ,pstar)
fprintf('   Shadow Price  %5.2f\n\n',lstar)

% Check model derivatives
dpcheck(model,sstar,qstar);

% Solve collocation equation
[v,q,c,sr,vr,qr,resid] = dpsolve(model,basis);

% Plot optimal policy
figure
hold on
plot(sr,qr)
plotbullet(sstar,qstar)
plotvdash(sstar,qstar)
plothdash(sstar,qstar)
plottext(sstar+0.02,[],'$s^*$')
plottext([],qstar,'$q^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Harvest Policy')
xlabel('Stock')
ylabel('Quantity Harvested')

% Plot value function
figure
hold on
plot(sr,vr)
plotvdash(sstar,vstar)
plothdash(sstar,vstar)
plotbullet(sstar,vstar)
plottext(sstar+0.01,[],'$s^*$')
plottext([],vstar,'$v^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Value Function')
xlabel('Stock')
ylabel('Social Welfare')

% Plot shadow price function
figure
hold on
pr = funeval(c,basis,sr,1);
plot(sr,pr)
plotvdash(sstar,lstar)
plothdash(sstar,lstar)
plotbullet(sstar,lstar)
plottext(sstar+0.02,[],'$s^*$')
plottext([],lstar,'$\lambda^*$')
xtickformat('%.1f')
ytickformat('%.2f')
title('Shadow Price Function')
xlabel('Stock')
ylabel('Shadow Price')

% Plot collocation equation residual
figure
hold on
plot(sr,resid)
plothdash([],0)
xtickformat('%.1f')
ytickformat('%.1f')
title('Bellman Equation Residual')
xlabel('Stock')
ylabel('Residual')


%% BASE CASE MODEL SIMULATION

% Initialize simulation
sinit = smin;

% Simulate model
[ssim,qsim] = dpsimul(model,basis,nper,sinit,[],sr,vr,qr);

% Plot simulated state and policy paths
figure
hold on
plot(0:nper-1,ssim,0:nper-1,qsim)
plothdash([],sstar,'b')
plothdash([],qstar,'r')
plottext([],sstar,'$s^*$')
plottext([],qstar,'$q^*$')
title('Simulated Stock and Harvest')
xlabel('Period')
ylabel('Quantity')
legend('Stock','Harvest')


%% PARAMETRIC SENSITIVITY ANALYSIS

fprintf('Varying Relative Risk Aversion\n')
gammamin = 0.0;
gammamax = 0.2;
gammaplot = nodeunif(nplot,gammamin,gammamax);
pstarplot = qstar.^(-gammaplot);
lstarplot = pstarplot-kappa;
figure
plot(gammaplot,lstarplot)
xticks(gammamin:0.1:gammamax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Relative Risk Aversion')
ylabel('Steady-State Shadow Price')

fprintf('Varying Discount Factor\n')
deltamin = 0.85;
deltamax = 0.95;
deltaplot = nodeunif(nplot,deltamin,deltamax);
sstarplot = (alpha^2-1./deltaplot.^2)/(2*beta);
qstarplot = sstarplot - (deltaplot*alpha-1)./(deltaplot*beta);
pstarplot = qstarplot.^(-gamma);
lstarplot = pstarplot-kappa;
figure
plot(deltaplot,sstarplot)
xticks(deltamin:0.05:deltamax)
xtickformat('%.2f')
ytickformat('%.2f')
xlabel('Discount Factor')
ylabel('Steady-State Stock')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,s,q,~,~,~,alpha,beta,gamma,kappa)
switch flag
  case 'b'      % bounds
    out1 = zeros(size(s));
    out2 = s;
    out3 = [];
  case 'f'      % reward
    out1 = (q.^(1-gamma))/(1-gamma)-kappa*q;
    out2 = q.^(-gamma)-kappa;
    out3 = -gamma*q.^(-gamma-1);
  case 'g'      % transition
    out1 = alpha*(s-q) - 0.5*beta*(s-q).^2;
    out2 = -alpha + beta*(s-q);
    out3 = -beta*ones(size(s));
end
end