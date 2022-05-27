%% DEMDP01 Timber Harvesting Model
%
% Profit maximizing a commercial tree stand owner must decide when to
% clear-cut and replant.

% States
%     s       stand biomass
% Actions
%     j       clear cut/replant (2), not clear cut (1)
% Parameters
%     price   unit price of biomass
%     kappa   clearcut-replant cost
%     sseed   biomass of seedling
%     scap    stand carrying capacity
%     delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'spli';                     % basis function type
n         = 200;                        % degree of approximation
nper      = 31;                         % number of periods simulated
time      = 0:nper-1;                   % periods simulated
nplot     = 21;                         % number of parameter values plotted

% Base case model parameters
price = 1.0;                           	% price of biomass
kappa = 0.2;                            % clearcut-replant cost
gamma = 0.05;                           % biomass growth parameter
scap  = 1.0;                            % stand carrying capacity
delta = 0.9;                            % discount factor

% Timber growth curve
h = @(s) s + gamma*(scap-s);
sseed = h(0);

% Refined set of state nodes
smax = 1.2*scap;
sr = nodeunif(100,0,smax);

% Plot timber growth curve
figure
hold on
plot(sr,h(sr))
plot(sr,sr,'k:')
plotbullet(scap,scap)
plotvdash(scap,scap)
plottext(0.15,0,'$45^\circ$','right','bottom',12)
plottext(smax,0,'$s$','center','top')
plottext(scap,0,'$\bar s$','center','top')
plottext(0,sseed,'$s_0$','right','middle',16)
plottext(0,smax,'$h(s)$','right','top',16)
axis square
xticks([])
yticks([])


%% BASE CASE MODEL SOLUTION - TWO COLLOCATION NODES

% Collocation nodes
snodes = [0.2;0.8]; 

% Refined set of state nodes
sr = nodeunif(100,0,scap);

% Value function approximants and residual
vhat  = @(c,s) c(1)+c(2)*s;
vhat0 = @(c,s)               delta*vhat(c,h(s));
vhat1 = @(c,s) price*s-kappa+delta*vhat(c,sseed);
resid = @(c,s) vhat(c,s) - max(vhat0(c,s),vhat1(c,s));

% Solve collocation equation with two collocation nodes
c = zeros(2,1);
c = broyden(resid,c,snodes);

% Critical biomass and value
scrit = interp1(vhat1(c,sr)-vhat0(c,sr),sr,0);
vcrit = vhat0(c,scrit);

% Compute optimal rotation cycle
s = sseed;
t = 1;
while s<scrit
  t = t+1;
  s = h(s);
end

% Output
fprintf('Two Collocation Nodes\n')
fprintf('Critical Biomass           %9.3f\n',scrit) 
fprintf('Mean Annual Harvest        %9.3f\n',s/t) 
fprintf('Rotation Cycle in Years    %9i\n\n',t) 

% Plot action-contingent value functions
figure
hold on
plot(sr,[vhat0(c,sr) vhat1(c,sr)])
plotvdash(scrit,vcrit)
plotbullet(scrit,vcrit)
plottext(scrit+0.01,-0.2,'$s^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Action-Contingent Value Functions - Two Collocation Nodes')
xlabel('Biomass')
ylabel('Value of Stand')
legend('Grow','Clear-Cut')

% Plot residual
figure
hold on
residpct = 100*resid(c,sr)./vhat(c,sr);
plot(sr,residpct)
plot(sr,0*sr,'k-','LineWidth',1)
plotvdash(snodes(1),0)
plotvdash(snodes(2),0)
plotvdash(scrit,max(residpct))
plottext(scrit+0.01,[],'$s^*$')
plottext(snodes(1)+0.01,[],'$s_1$')
plottext(snodes(2)+0.01,[],'$s_2$')
xtickformat('%.1f')
ytickformat('%.0f%%')
title('Bellman Equation Residual - Two Collocation Nodes')
xlabel('Biomass')
ylabel('Percent Residual')


%% BASE CASE MODEL SOLUTION

% Approximation structure
smax = 1.2*scap;
[basis,~,snodes] = fundefn(basistype,n,0,smax);

% Value function approximants and residual
vhat  = @(c,s) funeval(c,basis,s);
vhat0 = @(c,s)               delta*vhat(c,h(s));
vhat1 = @(c,s) price*s-kappa+delta*vhat(c,sseed);
resid = @(c,s) vhat(c,s) - max(vhat0(c,s),vhat1(c,s));

% Solve collocation equation
c = zeros(n,1);
c = broyden(resid,c,snodes);

% Critical biomass and value
scrit = interp1(vhat1(c,sr)-vhat0(c,sr),sr,0);
vcrit = vhat0(c,scrit);

% Compute optimal rotation cycle
s = sseed;
t = 1;
while s<scrit
  t = t+1;
  s = h(s);
end

% Output
fprintf('Many Collocation Nodes\n')
fprintf('Critical Biomass           %9.3f\n',scrit) 
fprintf('Mean Annual Harvest        %9.3f\n',s/t) 
fprintf('Rotation Cycle in Years    %9i\n\n',t) 

% Plot action-contingent value functions
figure
hold on
plot(sr,[vhat0(c,sr) vhat1(c,sr)])
plotvdash(scrit,vcrit)
plotbullet(scrit,vcrit)
plottext(scrit+0.01,[],'$s^*$')
xtickformat('%.1f')
ytickformat('%.1f')
title('Action-Contingent Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')
legend('Grow','Clear-Cut')

% Plot residual
figure
hold on
residpct = 100*resid(c,sr)./max(vhat0(c,sr),vhat1(c,sr));
plot(sr,residpct)
plothdash([],0)
plotvdash(scrit,max(residpct))
plottext(scrit+0.01,[],'$s^*$')
xtickformat('%.1f')
ytickformat('%.2f%%')
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')


%% BASE CASE MODEL SIMULATION

% Initialize simulation
s = sseed;                                  % initial biomass

% Simulate model
ssim = zeros(nper,1);
for ip=1:nper
  ssim(ip) = s;
  if s<scrit
    s = h(s);
  else
    s = sseed;
  end
end

% Plot state path
figure
hold on
bar(time,ssim)
plothdash([],scrit)
plottext(0,scrit,'$s^*$','left','bottom')
ytickformat('%.1f')
title('Simulated Biomass')
xlabel('Period')
ylabel('Biomass')


%% PARAMETRIC SENSITIVITY ANALYSIS

% Intialization
scritplot = zeros(nplot,1);

fprintf('Varying Price of Biomass\n')
priceplot = nodeunif(nplot,0.6,1.4);
c = zeros(n,1);
for ip=1:nplot
  vhat0 = @(c,s)                       delta*vhat(c,h(s));
  vhat1 = @(c,s) priceplot(ip)*s-kappa+delta*vhat(c,sseed);
  resid = @(c,s) vhat(c,s) - max(vhat0(c,s),vhat1(c,s));
  c = broyden(resid,c,snodes);
  scritplot(ip) = interp1(vhat1(c,sr)-vhat0(c,sr),sr,0);
end
figure
plot(priceplot,scritplot)
xtickformat('%.1f')
ytickformat('%.2f')
xlabel('Price of Biomass')
ylabel('Critical Biomass')

fprintf('Varying Clearcut-Replant Cost\n')
kappaplot = nodeunif(nplot,0.0,0.4);
c = zeros(n,1);
for ip=1:nplot
  vhat0 = @(c,s)                       delta*vhat(c,h(s));
  vhat1 = @(c,s) price*s-kappaplot(ip)+delta*vhat(c,sseed);
  resid = @(c,s) vhat(c,s) - max(vhat0(c,s),vhat1(c,s));
  c = broyden(resid,c,snodes);
  scritplot(ip) = interp1(vhat1(c,sr)-vhat0(c,sr),sr,0);
end
figure
plot(kappaplot,scritplot)
xtickformat('%.1f')
ytickformat('%.1f')
xlabel('Clearcut-Replant Cost')
ylabel('Critical Biomass')


%% SAVE FIGURES
printfigures(mfilename)