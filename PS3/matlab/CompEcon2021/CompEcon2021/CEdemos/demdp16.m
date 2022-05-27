%% DEMDP16 Savings and Insurance
%
% An infinitely-lived agent begins each period endowed with predetermined
% wealth w, which he must allocate among consumption, savings, and
% purchases of insurance. 

% States
%     w       stock of wealth
% Actions
%     s       savings
%     x       insurance coverage
% Parameters
%     alpha   relative risk aversion
%     rho     discount rate
%     r       interest rate on savings
%     ybar    income mean
%     L       magnitude of loss
%     p       probability of loss
%     load    premium load

% Preliminary tasks
deminit(mfilename)



%% FORMULATION

% Numerical control parameters
n     = 200;                            % degree of approximation
nper  = 1000;                           % number of periods simulated
xlm   = [0 2.5];                        % x-axis plot limits
ylm   = [0 0.7];                        % y-axis plot limits
nplot = 15;                             % number of parameter values plotted

% Base case model parameters
smax  = 3.5;                            % maximum savings - 0 or 2.5-3.5
xmax  = 1;                              % maximum coverage - 0 or 1
alpha = 3.0;                           	% relative risk aversion
rho   = 0.10;                           % discount rate
r     = 0.08;                           % interest rate on savings
ybar  = 1.0;                            % income mean
L     = 0.5;                            % magnitude of insurable loss
p     = 0.2;                            % probability of insurable loss
load  = 0.2;                            % premium load
if r>rho, warning('Interest rate exceeds discount rate'), end

% Derived parameters
delta = 1/(1+rho);                      % discount factor
prem = (1+load)*delta*p*L;              % premium rate

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {alpha ybar prem r smax xmax};% function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 2;                           % dimension of continuous action
model.e  = [0;L];                       % shocks
model.w  = [1-p;p];                     % shock probabilities

% Approximation structure
wmin  = ybar-L;                         % minimum state
wmax  = ybar+(1+r)*smax;                % maximum state
basis = fundefn('spli',n,wmin,wmax);    % basis functions
if wmin<0,error('Negative wealth in second period'), end

% Check model derivatives
dpcheck(model,1,[0 0])


%% BASE CASE MODEL SOLUTION

% Both savings and insurance available
[~,wA,vA,xA,residA,sseA] = solve(ybar,L,p,model,basis,nper,r,prem);

% Only insurance available, no savings
model.params = {alpha ybar prem r 0 xmax};
[~,wI,vI,xI,residI,sseI] = solve(ybar,L,p,model,basis,nper,r,prem);

% Only savings available, no insurance
model.params = {alpha ybar prem r smax 0};
[~,wS,vS,xS,residS,sseS] = solve(ybar,L,p,model,basis,nper,r,prem);

% Neither savings nor insurance available
model.params = {alpha ybar prem r 0 0};
[~,w0,v0,~,resid0,~] = solve(ybar,L,p,model,basis,nper,r,prem);

% Print ergodic means of savings and insurance expenditures
fprintf('\n')
fprintf('Ergodic Mean Expenditures as Proportion of Wealth\n')
fprintf('               Savings and  Insurance    Savings\n')
fprintf('                Insurance     Only        Only\n')
fprintf('%-10s%8.3f%11.3f%12.3f\n','Savings        ',sseA(1),sseI(1),sseS(1))
fprintf('%-10s%8.3f%11.3f%12.3f\n\n','Insurance      ',sseA(2),sseI(2),sseS(2))

% Plot optimal savings and insurance policy
figure
hold on
plot(wA,xA(:,1)./wA,'b')
plot(wA,prem*xA(:,2)./wA,'r')
plot(wS,xS(:,1)./wS,'b:')
plot(wI,prem*xI(:,2)./wI,'r:')
xlim([0 3])
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Savings and Insurance Expenditures')
xlabel('Wealth')
ylabel('Expenditure as Proportion of Wealth')
legend('Savings', 'Insurance', ...
       'Savings w/o Insurance','Insurance w/o Savings')

% Plot optimal insurance coverage policy
figure
hold on
plot(wA,xA(:,2),'r')
plot(wI,xI(:,2),'r:')
xlim([0 3])
ylim([0 1])
xtickformat('%.1f')
ytickformat('%.1f')
title('Optimal Insurance Coverage')
xlabel('Wealth')
ylabel('Coverage')
legend('With Savings', 'Without Savings')

% Plot value function
figure
plot(wA,vA,wS,vS,wI,vI,w0,v0)
title('Value Function')
xlabel('Wealth')
ylabel('Lifetime Utility')
legend('Savings and Insurance','Savings Only','Insurance Only','Neither')

% Plot residuals
figure
plot(wA,residA,wS,residS,wI,residI,w0,resid0)
ytickformat('%.1f')
title('Residual')
xlabel('Wealth')
ylabel('Residual')
legend('Savings and Insurance','Savings Only','Insurance Only','Neither')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialize
covplot = zeros(1,nplot);
sseplot = zeros(2,nplot);

fprintf('Varying Premium Load\n')
modelplot = model;
loadmin = -0.5;
loadmax =  0.5;
loadplot = nodeunif(nplot,loadmin,loadmax);
for ip=1:nplot
  premplot = (1+loadplot(ip))*delta*p*L;
  modelplot.params = {alpha ybar premplot r smax xmax};
  [~,~,~,~,~,sseplot(:,ip),covplot(ip)] = solve(ybar,L,p,modelplot,basis,nper,r,premplot);
end
figure
plot(100*loadplot,sseplot)
xtickformat('%.0f%%')
ytickformat('%.1f')
title('Ergodic Mean Savings and Insurance Expenditures vs. Premium Load')
xlabel('Premium Load')
ylabel('Proportion of Wealth')
legend('Savings','Insurance Expenditure')
figure
plot(100*loadplot,covplot)
xtickformat('%.0f%%')
ytickformat('%.1f')
title('Ergodic Mean Coverage vs. Premium Load')
xlabel('Premium Load')
ylabel('Proportion of Loss Insured')

fprintf('Varying Interest Rate\n')
modelplot = model;
rmin = 0.00;
rmax = 0.08;
rplot = nodeunif(nplot,rmin,rmax);
for ip=1:nplot
  modelplot.params = {alpha ybar prem rplot(ip) smax xmax};
  [~,~,~,~,~,sseplot(:,ip),covplot(ip)] = solve(ybar,L,p,modelplot,basis,nper,rplot(ip),prem);
end
figure
plot(100*rplot,sseplot)
xtickformat('%.0f%%')
ytickformat('%.2f')
title('Ergodic Mean Savings and Insurance Expenditures vs. Interest Rate')
xlabel('Interest Rate')
ylabel('Proportion of Wealth')
legend('Savings','Insurance Expenditure')
figure
plot(100*rplot,covplot)
xtickformat('%.0f%%')
ytickformat('%.1f')
title('Ergodic Mean Coverage vs. Interest Rate')
xlabel('Interest Rate')
ylabel('Proportion of Loss Insured')


%% SAVE FIGURES
printfigures(mfilename)


%% ANCILLARY FUNCTION
%
% Solves and simulates model

function [c,wr,vr,xr,resid,sse,cov] = solve(ybar,L,p,model,basis,nper,r,prem)
% Solve collocation equation
optset('dpsolve','output',0)
optset('dpsolve','maxit',15)
optset('dpsolve','algorithm','funcit')
[v,x] = dpsolve(model,basis);
optset('dpsolve','algorithm','newton')
[~,~,c,wr,vr,xr,resid] = dpsolve(model,basis,v,x);
% Simulate model
if nper>0
  rng('default')
  wsim  = zeros(nper,1);
  xsim  = zeros(nper,2);
  i = rand(nper,1)>1-p;
  wtmp  = 1;
  for ip=1:nper
    wsim(ip) = wtmp;
    xsim(ip,:) = interp1(wr,xr,wtmp);
    wtmp = ybar + (1+r)*xsim(ip,1) - (1-xsim(ip,2))*i(ip)*L;
  end
  sse = [mean(xsim(:,1)./wsim); mean(prem*xsim(:,2)./wsim)];
  cov = mean(xsim(:,2));
end
end


%% DPSOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,w,x,~,~,l,alpha,ybar,prem,r,smax,xmax)
n  = length(w);
switch flag
  case 'b'      % bounds
    out1 = zeros(n,2);
    out2 = [smax+zeros(n,1) xmax+zeros(n,1)];
    out3 = [];
  case 'f'      % reward
    out2 = zeros(n,2);
    out3 = zeros(n,2,2);
    c = w - x(:,1) - prem*x(:,2);
    dudc1 = c.^(-alpha);
    dudc2 = -alpha*c.^(-alpha-1);
    dcdx1 = -1;
    dcdx2 = -prem;
    out1        = (c.^(1-alpha))/(1-alpha);
    out2(:,1)   = dudc1*dcdx1;
    out2(:,2)   = dudc1*dcdx2;
    out3(:,1,1) = dudc2*dcdx1*dcdx1;
    out3(:,1,2) = dudc2*dcdx1*dcdx2;
    out3(:,2,1) = dudc2*dcdx2*dcdx1;
    out3(:,2,2) = dudc2*dcdx2*dcdx2;
  case 'g'      % transition
    out2 = zeros(n,1,2);
    out3 = zeros(n,1,2,2);
    out1 = ybar + (1+r)*x(:,1) - (1-x(:,2)).*l;
    out2(:,1,1) = (1+r)*ones(n,1);
    out2(:,1,2) = l;
end
end