%% DEMDP15 Lifecycle Consumption-Saving Model
%
% In the lifecycle consumption-saving model, an agent works from period 1
% to period T, earning a known income y(t) in each period t, a portion g of
% which he is required to pay into a pension fund. In this code, the Euler
% conditions are solved as a nonlinear complementarity problem using
% ncpsolve.

% Endogenous Variable
%     x       net ending liquid assets (x<0 implies borrowing)
% Parameters
%     g       pension contribution rate (proportion of income)
%     r       interest rate earned by assets and paid on debt
%     alpha   agent's constant relative risk aversion
%     T       number of working periods
%     tmax    period of maximum income
%     ymax    maximum income (income at t=1 is 1)
%     k       borrowing limit as a proportion of income
%     delta   agent's subjective discount factor
% Derived Parameters
%     R       gross interest rate (1+r)
%     y       income per period
%     d       disposable income per period

% Preliminary tasks
close all
clear all


%% FORMULATION

% Numerical control parameters
nplot = 25;         % number of parameter values plotted

% Base case model parameters
g     = 0.0;        % pension contribution rate (proportion of income)
r     = 0.1;        % interest rate earned by assets, paid on debt
alpha =   4;        % agent's relative risk aversion
T     =  40;        % number of periods until retirement
k     = 0.5;        % borrowing limit as a proportion of income
delta = 0.9;        % agent's subjective discount factor
tmax  = 0.75*T;     % period of maximum income
ymax  = 1.5;        % maximum income (income=1 at t=0)

% Derived parameters
t = (0:T)';         % time periods
R = 1+r;            % gross interest rate

% Utility function and derivative
u    = @(c) c.^(1-alpha)/(1-alpha);
uder = @(c) c.^(-alpha);

% Income stream
c = 2*(ymax-1)/tmax^2;
b = tmax*c;
y = 1 + b*t - 0.5*c*t.^2;

% Plot income stream
figure
plot(t,y)
ytickformat('%.1f')
title('Agent''s Income Stream')
xlabel('Period')
ylabel('Income')

% Disposable income stream
d = (1-g)*y;

% Value of pension at retirement
G = g*sum(y.*(R.^(T-t+1)));


%% BASE CASE MODEL SOLUTION

% Solve Euler conditions as nonlinear complementarity problem
[x,e] = ncpsolve(@F,-k*d,inf,0.2*d,d,T,delta,R,r,G,uder);

% Consumption
c = R*[0;x(1:T)] + d - x;
if any(c<0), disp('WARNING: Negative Consumption'), end

% Plot conumption and liquid assets over time
figure
hold on
plot(t,[c x])
plothdash([],0)
title('Consumption and Liquid Assets')
xlabel('Period')
ylabel('Amount')
legend('Consumption','Liquid Assets')

% Compute equivalent constant lifetime consumption 
U = sum((delta.^t).*u(c)) + (delta^T/(1-delta))*u(r*(G+R*x(T+1)));
clife = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));

% Print output
fprintf('\n')
fprintf('Liquid Assets at Retirement               %5.2f\n',G+R*x(T))
fprintf('Equivalent Lifetime Constant Consumption  %5.2f\n',clife)
fprintf('Norm of Euler Equation                 %5.2e\n\n',norm(e))


%% PARAMETRIC SENSITIVITY ANALYSIS - BORROWING LIMITS

% Solve for Different Borrowing Limits
kvec = [0 1 2];
N = length(kvec);

% Solve Euler Conditions
x = zeros(T+1,N);
for i=1:N
  x(:,i) = ncpsolve(@F,-kvec(i)*d,inf,0.1*d,y,T,delta,R,r,G,uder);
end
s = R*[zeros(1,N);x(1:T,:)];

% Plot liquid assets
figure
hold on
plot(t,s)
plothdash([],0)
title('Liquid Assets - Different Borrowing Limits')
xlabel('Period')
ylabel('Liquid Assets')
legend('$k=0$','$k=1$','$k=2$')


%% PARAMETRIC SENSITIVITY ANALYSIS - INTEREST RATES

% Solve for Different Interest Rates
rvec = [0.05 0.10 0.15];
Rvec = 1+rvec;
Gvec = g*sum(y.*(Rvec.^(T-t+1)));
N = length(rvec);

% Solve Euler Conditions
x = zeros(T+1,N);
for i=1:N
  x(:,i) = ncpsolve(@F,-k*d,inf,0.1*d,y,T,delta,Rvec(i),rvec(i),Gvec(i),uder);
end
s = Rvec.*[zeros(1,N);x(1:T,:)];

% Plot liquid assets
figure
hold on
plot(t,s)
plothdash([],0)
title('Liquid Assets - Different Interest Rates')
xlabel('Period')
ylabel('Liquid Assets')
legend('$r=0.05$','$r=0.10$','$r=0.15$')


%% PARAMETRIC SENSITIVITY ANALYSIS - PENSION CONTRIBUTION RATES

% Solve for Different Pension Contribution Rates
gvec = [0.0 0.1 0.2];
Gvec = gvec*sum(y.*(R.^(T-t+1)));
N = length(gvec);

% Solve Euler Conditions
x = zeros(T+1,N);
for i=1:N
  x(:,i) = ncpsolve(@F,-k*d,inf,0.1*d,y,T,delta,R,r,Gvec(i),uder);
end
s = R*[zeros(1,N);x(1:T,:)];

% Plot liquid assets
figure
hold on
plot(t,s)
plothdash([],0)
title('Liquid Assets - Different Pension Fund Contribution Rates')
xlabel('Period')
ylabel('Liquid Assets')
legend('$g=0.0$','$g=0.1$','$g=0.2$')


%% OTHER PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
clife = zeros(nplot,1);

fprintf('Varying Interest Rate\n')
rmin = 0.01;
rmax = 0.30;
rplot = nodeunif(nplot,rmin,rmax);
for ip=1:nplot
  Rplot = 1+rplot(ip);
  Gplot = g*sum(y.*(Rplot.^(T-t+1)));
  x = ncpsolve(@F,-k*d,inf,0.2*d,d,T,delta,Rplot,rplot(ip),Gplot,uder);
  c = Rplot*[0;x(1:T)] + d - x;
  U = sum((delta.^t).*u(c)) + (delta^T/(1-delta))*u(rplot(ip)*(Gplot+Rplot*x(T+1)));
  clife(ip) = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));
end
figure
plot(rplot,clife)
xlim([0.0 0.3])
ylim([0.8 1.2])
xticks(0.0:0.1:0.3)
yticks(0.8:0.1:1.2)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Interest Rate')
ylabel('Equivalent Constant Consumption')

fprintf('Varying Pension Contribution Rate\n')
gmin = 0.0;
gmax = 0.3;
gplot = nodeunif(nplot,gmin,gmax);
for ip=1:nplot
  Gplot = gplot(ip)*sum(y.*(R.^(T-t+1)));
  dplot = (1-gplot(ip))*y;
  x = ncpsolve(@F,-k*dplot,inf,0.2*dplot,dplot,T,delta,R,r,Gplot,uder);
  c = R*[0;x(1:T)] + dplot - x;
  U = sum((delta.^t).*u(c)) + (delta^T/(1-delta))*u(r*(Gplot+R*x(T+1)));
  clife(ip) = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));
end
figure
plot(gplot,clife)
ylim([0.8 1.2])
xticks(0.0:0.1:0.3)
yticks(0.8:0.1:1.2)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Pension Contribution Rate')
ylabel('Equivalent Constant Consumption')

fprintf('Varying Borrowing Limit\n')
kmin = 0.0;
kmax = 1.0;
kplot = nodeunif(nplot,kmin,kmax);
for ip=1:nplot
  x = ncpsolve(@F,-kplot(ip)*d,inf,0.2*d,d,T,delta,R,r,G,uder);
  c = R*[0;x(1:T)] + d - x;
  U = sum((delta.^(t-1)).*u(c)) + (delta^T/(1-delta))*u(r*(G+R*x(T+1)));
  clife(ip) = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));
end
figure
plot(kplot,clife)
ylim([0.8 1.2])
yticks(0.8:0.1:1.2)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Borrowing Limit')
ylabel('Equivalent Constant Consumption')

fprintf('Varying Relative Risk Aversion\n')
alphamin = 1;
alphamax = 5;
alphaplot = nodeunif(nplot,alphamin,alphamax);
for ip=1:nplot
  uplot    = @(c) c.^(1-alphaplot(ip))/(1-alphaplot(ip));
  uderplot = @(c) c.^(-alphaplot(ip));
  x = ncpsolve(@F,-k*d,inf,0.5*d,d,T,delta,R,r,G,uderplot);
  c = R*[0;x(1:T)] + d - x;
  U = sum((delta.^(t-1)).*uplot(c)) + (delta^T/(1-delta))*uplot(r*(G+R*x(T+1)));
  clife(ip) = ((1-delta)*(1-alphaplot(ip))*U)^(1/(1-alphaplot(ip)));
end
figure
plot(alphaplot,clife)
ylim([0.8 1.2])
xticks(1:1:5)
yticks(0.8:0.1:1.2)
ytickformat('%.1f')
xlabel('Relative Risk Aversion')
ylabel('Equivalent Constant Consumption')


%% SAVE FIGURES
printfigures(mfilename)


%% FUNCTION FILE
function [f,J] = F(x,y,T,delta,R,r,G,uder)
ud = [uder(R*[0;x(1:T)]+y-x); (r/(1-delta))*uder(r*(G+R*x(T+1))) ];
t = (1:T+1)';
f = -ud(t) + delta*R*ud(t+1);
J = [];
end