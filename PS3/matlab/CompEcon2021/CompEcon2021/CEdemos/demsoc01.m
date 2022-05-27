%% DEMSOC01 Stochastic Consumption-Investment Model
%
% Agent must decide the rate at which to consume and invest in risky asset.

% State
%     w       wealth
% Control
%     q       consumption rate
% Parameters
%     r       expected rate of return on asset
%     theta   relative risk aversion
%     sigma   volatility or return on asset
%     rho     continuous discount rate

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
winit = 1;              % initial wealth
T     = 1;              % time horizon
N     = 1000;           % number of time intervals
nrep  = 10000;          % number of replications

% Base case model parameters
r     = 0.08;           % expected rate of return on asset
theta = 2.0;            % relative risk aversion
sigma = 0.05;           % volatility or return on asset
rho   = 0.05;           % continuous discount rate


%% BASE CASE MODEL SOLUTION

A = (rho-r*(1-theta)+0.5*theta*(1-theta)*sigma^2)/theta;
if A<=0
  fprinft('\nInviable Parameters\n')
end


%% BASE CASE MODEL SIMULATION

% Time discretization
h = T/N;                % length of time intervals
t = (0:h:T)';           % time nodes

% Preallocate output array
wsim = zeros(N+1,nrep);

% Simulate model
z = sqrt(h)*randn(N,nrep)*sigma;
w = winit*ones(1,nrep);
for i=0:N
  wsim(i+1,:,:) = w;
  if i<N
    w = w + (r-A)*w*h + w.*z(i+1,:);
  end
end

% Plot optimal wealth paths
figure
hold on
plot(t,wsim(:,1:3),t,mean(wsim,2),'k')
xtickformat('%.1f')
ytickformat('%.2f')
title('Simulated and Expected Wealth')
xlabel('Time')
ylabel('Wealth')


%% SAVE FIGURES
printfigures(mfilename)