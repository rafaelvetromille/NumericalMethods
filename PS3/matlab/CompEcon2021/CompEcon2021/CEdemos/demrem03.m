%% DEMREM03 Government Price Support Model
%
% A government buffer stock authority defends a minimum market price
% through the acquisition of excess supply.

% States
%     s         supply
% Responses
%     a         acreage planted
%     z         government stocks
% Parameters
%     zmax      maximum stocks
%     pbar      government support price
%     gamma     inverse consumption demand elasticity
%     beta      acreage supply elasticity
%     sigma     yield volatility
%     delta     discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
n     = 501;                                        % degree of approximation
smin  = 0.5;                                        % minimum supply
smax  = 2.5;                                        % maximum supply
m     = 11;                                         % number of yield shocks
nper  = 30;                                         % number of periods simulated
nrep  = 10000;                                      % number of replications
nplot = 51;                                         % number of parameter values plotted

% Base case model parameters
zmax  = 0.5;                                        % maximum government stocks
pbar  = 1.0;                                        % government support price
gamma = 0.5;                                        % inverse consumption demand elasticity
beta  = 1.0;                                        % acreage supply elasticity
sigma = 0.2;                                        % yield volatility
delta = 0.9;                                        % discount factor

% Discretized yield distribution
% [y,w] = qnwlogn(m,1,sigma^2);
[y,w] = qnwlognunif(m,1,sigma^2);

% Approximation structure
s = nodeunif(n,smin,smax);    


%% BASE CASE MODEL SOLUTION

% Critical supply levels
scrit = pbar^(-1/gamma);
p   = @(s) s.^(-gamma);
zeq = @(s) min(max(s-scrit,0),zmax);
peq = @(s) p(s-zeq(s));
z = zeq(s);

% Solve rational expectations equilibrium acreage using broyden
a = broyden(@f,ones(n,1),z,m,n,y,w,delta,beta,peq);
% a = ones(n,1);

% Plot equilibrium planged acreage
figure
hold on
plot(s,a)
plotvdash(scrit,max(a))
plotvdash(scrit+zmax,min(a))
plottext(scrit+0.02,[],'$s_1^*$')
plottext(scrit+zmax+0.02,[],'$s_2^*$')
xtickformat('%.1f')
ytickformat('%.2f')
title('Equilibrium Planted Acreage')
xlabel('Supply')
ylabel('Acreage')

% Plot equilibrium government stocks
figure
hold on
plot(s,z)
plotvdash(scrit,0)
plotvdash(scrit+zmax,zmax)
plottext(scrit+0.02,[],'$s_1^*$')
plottext(scrit+zmax+0.02,[],'$s_2^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Equilibrium Government Stocks')
xlabel('Supply')
ylabel('Stocks')

% Plot equilibrium market price function
figure
hold on
plot(s,[peq(s) p(s)])
plot([0 smax],[pbar pbar],'k:')
ylim([0.5 1.5])  
plotvdash(scrit,pbar)
plotvdash(scrit+zmax,pbar)
plottext(scrit+0.02,[],'$s_1^*$')
plottext(scrit+zmax+0.02,[],'$s_2^*$')
ytickformat('%.1f')
title('Equilibrium Market Price')
xlabel('Supply')
ylabel('Price')
legend('Total Demand','Consumption Demand')


%% BASE CASE MODEL SIMULATION

% Generate random yields
rng('default')
ysim = randlogn(1,sigma^2,nrep,nper+1);

% Simulate model
[ssim,zsim,psim,rsim,asim] = simul(nrep,nper,zmax,ysim,gamma,pbar,s,a);

% Plot simulated and expected supply
figure
plot(0:nper,ssim(1:3,:),0:nper,mean(ssim),'k')
ytickformat('%.1f')
title('Simulated and Expected Supply')
xlabel('Period')
ylabel('Supply')

% Plot simulated and expected planted acreage
figure
plot(0:nper,asim(1:3,:),0:nper,mean(asim),'k')
title('Simulated and Expected Planted Acreage')
ytickformat('%.2f')
xlabel('Period')
ylabel('Acreage')

% Plot simulated and expected government stocks
figure
plot(0:nper,zsim(1:3,:),0:nper,mean(zsim),'k')
ytickformat('%.1f')
title('Simulated and Expected Government Stocks')
xlabel('Period')
ylabel('Stocks')

% Plot simulated and expected market price
figure
plot(0:nper,psim(1:3,:),0:nper,mean(psim),'k')
ylim([0.5 1.5])
ytickformat('%.1f')
title('Simulated and Expected Market Price')
xlabel('Period')
ylabel('Price')

% Ergodic moments
savg = mean(ssim(:));
aavg = mean(asim(:));
zavg = mean(zsim(:));
pavg = mean(psim(:));
ravg = mean(rsim(:));
sstd = std(ssim(:));
astd = std(asim(:));
zstd = std(zsim(:));
pstd = std(psim(:));
rstd = std(rsim(:));

% Print ergodic moments
fprintf('\n')
fprintf('                 Ergodic      Ergodic\n') 
fprintf('                  Mean     Std Deviation\n') 
fprintf('Supply           %5.3f         %5.3f\n',  [savg sstd])
fprintf('Acreage Planted  %5.3f         %5.3f\n',  [aavg astd])
fprintf('Ending Stocks    %5.3f         %5.3f\n',  [zavg zstd])
fprintf('Market Price     %5.3f         %5.3f\n',  [pavg pstd])
fprintf('Producer Revenue %5.3f         %5.3f\n\n',[ravg rstd])

% Plot ergodic supply distribution
[qq,ss] = ksdensity(ssim(:),'support','positive','bandwidth',0.05);
figure
plot(ss,qq)
xlim([0 2])  
xtickformat('%.1f')
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Density')

% Plot ergodic supply distribution
figure
histogram(ssim(:),'Normalization','pdf')
xlim([0 2])  
xtickformat('%.1f')
title('Ergodic Supply Distribution')
xlabel('Supply')
ylabel('Density')

% Plot ergodic market price distribution
figure
histogram(psim(:),'Normalization','pdf')
xlim([0 3])  
xtickformat('%.1f')
title('Ergodic Market Price Distribution')
xlabel('Market Price')
ylabel('Density')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initialization
pavgplot = zeros(nplot,1);
pstdplot = zeros(nplot,1);
ravgplot = zeros(nplot,1);
rstdplot = zeros(nplot,1);

fprintf('Varying Government Support Price\n')
pbarmin = 0.0;
pbarmax = 1.5;
pbarplot = nodeunif(nplot,pbarmin,pbarmax);
for ip=1:nplot
  scrit = pbarplot(ip)^(-1/gamma);
  x = min(max(s-scrit,0),zmax);
  [ssim,xsim,psim,rsim] = simul(nrep,nper,zmax,ysim,gamma,pbarplot(ip),s,a);
  pavgplot(ip) = mean(psim(:,nper));
  pstdplot(ip) = std(psim(:,nper));
  ravgplot(ip) = mean(rsim(:,nper));
  rstdplot(ip) = std(rsim(:,nper));
end
figure
plot(pbarplot,pavgplot)
xlim([pbarmin pbarmax]) 
xticks(pbarmin:0.5:pbarmax)
xtickformat('%.1f')
ytickformat('%.3f')
xlabel('Government Support Price')
ylabel('Ergodic Mean Market Price')
figure
plot(pbarplot,pstdplot)
xlim([pbarmin pbarmax]) 
xticks(pbarmin:0.5:pbarmax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Government Support Price')
ylabel('Ergodic Standard Deviation of Market Price')
figure
plot(pbarplot,rstdplot)
xlim([pbarmin pbarmax]) 
xticks(pbarmin:0.5:pbarmax)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Government Support Price')
ylabel('Ergodic Standard Deviation of Revenue')


%% SAVE FIGURES
printfigures(mfilename)


%% FUNCTIONS

function fval = f(a,z,m,n,y,w,delta,beta,peq)
  sn   = z(:,ones(1,m))+a(:,ones(1,m)).*y(:,ones(1,n))';
  fval = a.^beta - delta*reshape(peq(sn(:)),n,m)*w;
end

function [ssim,zsim,psim,rsim,asim] = simul(nrep,nper,zmax,ysim,gamma,pbar,s,a)
ssim = zeros(nrep,nper+1);
zsim = zeros(nrep,nper+1);
asim = zeros(nrep,nper+1);
zz = zeros(nrep,1);
aa = ones(nrep,1);
scrit = pbar^(-1/gamma);
for ip=1:nper+1
  ss = zz+aa.*ysim(:,ip);
  zz = min(max(ss-scrit,0),zmax);
  aa = minterp(s,a,ss);
  asim(:,ip,:) = aa;
  ssim(:,ip,:) = ss;
  zsim(:,ip,:) = zz;
end
psim = (ssim-zsim).^(-gamma);
rsim = psim.*asim.*ysim;
end