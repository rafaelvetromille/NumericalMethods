%% DEMGAME02 Income Redistribution Game Model
%
% Two infinitely-lived agents who agree to share income risk must make 
% consumption and investment decisions. 

%    States
%        s       wealth, post income transfer
%    Actions
%        x       investment
%    Parameters
%        alpha   utility parameter
%        beta    production elasticity
%        gamma   capital survival rate
%        psi     risk sharing rate
%        sigma   production shock volatility
%        delta   discount factor

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Numerical control parameters
basistype = 'cheb';                     % basis function type
n    = [9 9];                           % degree of approximation
smin = [ 2  2];                         % minimum state
smax = [12 12];                         % maximum state
m    = [3 3];                           % number of shocks
nper = 30;                              % number of periods simulated
nrep = 5000;                            % number of replications

% Base case model parameters
alpha = [0.3 0.3];                      % utility parameter
beta  = [0.5 0.5];                      % production elasticity
gamma = [0.9 0.9];                      % capital survival rate
mu    = [0.9 1.1];                      % expected income
psi   = 0.1;                            % risk sharing rate
sigma = [0.2 0.2];                      % production shock volatility
evar  = diag(sigma.^2);                 % production shock variance matrix
delta = 0.9;                            % discount factor

% Shock distribution
[e,w] = qnwlogn(m,mu,evar);

% Model structure
clear model
model.func = @func;                     % model functions
model.params = {alpha beta gamma psi};  % other parameters
model.discount = delta;                 % discount factor
model.dp = 2;                           % number of players
model.ds = 2;                           % dimension of state variable s
model.dx = 1;                           % dimension of individual action variable x_i
model.e = e;                            % shocks
model.w = w;                            % probabilities

% Compute nonstochastic steady-state
xstar = ((1/delta-gamma)./beta).^(1./(beta-1));
sstar = gamma.*xstar+xstar.^beta;
ratio = xstar./sstar;
vstar = (1/(1-delta))*((sstar-xstar).^(1-alpha))./(1-alpha);
pstar = (sstar-xstar).^(-alpha);

% Approximation Structure
[basis,~,snodes] = fundefn(basistype,n,smin,smax);

% Initialize arrays
xinit = [snodes(:,1)*ratio(1) snodes(:,2)*ratio(2)];
vinit = [vstar(1)+pstar(1)*(snodes(:,1)-sstar(1)) vstar(2)+pstar(2)*(snodes(:,2)-sstar(2))] ;

% Check Model Derivatives
gamecheck(model,(smax+smin)/2,ratio.*(smax+smin)/2);


%% BASE CASE MODEL SOLUTION

% Solve Bellman Equation
optset('gamesolve','nres',3)
[v,x,c,sr,vr,xr,resid] = gamesolve(model,basis,vinit,xinit);

% Reshape Output for Plotting
nr  = n*3;
sr1 = reshape(sr(:,1),nr);
sr2 = reshape(sr(:,2),nr);
vr  = reshape(vr,[nr 2]);
xr  = reshape(xr,[nr 2]);
resid = reshape(resid,[nr 2]);

% Compute Shadow Prices
pr = zeros(prod(nr),2);
pr(:,1) = funeval(c(:,1),basis,sr,[1 0]);
pr(:,2) = funeval(c(:,2),basis,sr,[0 1]);
pr = reshape(pr,[nr 2]);

% Plot Optimal Policy
for ip=1:2
  figure
  mesh(sr1,sr2,xr(:,:,ip))
  xlim([0 15])
  ylim([0 15])
  title(['Optimal Investment Policy: Agent ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Investment')
end

% Plot Own Shadow Price Function
for ip=1:2
  figure
  mesh(sr1,sr2,pr(:,:,ip))
  xlim([0 15])
  ylim([0 15])
  ztickformat('%.1f')
  title(['Own Shadow Price: Agent ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Price')
end

% Plot Residual
for ip=1:2
  figure
  mesh(sr1,sr2,resid(:,:,ip))
  xlim([0 15])
  ylim([0 15])
  title(['Approximation Residual: Agent ' num2str(ip)])
  xlabel('Wealth 1')
  ylabel('Wealth 2')
  zlabel('Residual')
end


%% BASE CASE MODEL SIMULATION

% Generate random shocks
rng('default')
esim = randlogn(mu,evar,nrep,nper+1);

% Initialize simulation
nr  = n*3;
sr1 = reshape(sr(:,1),nr);
sr2 = reshape(sr(:,2),nr);
ss    = cell(2,1);
ss{1} = sr1(:,1);
ss{2} = sr2(1,:)';
sinit = 7.49*ones(nrep,2);

% Simulate Model
[ssim,xsim] = gamesimul(model,sinit,nper,ss,xr,esim);

% Plot simulated and expected wealth
for ip=1:2
  figure
  hold on
  plot(0:nper,ssim(1:3,:,ip))
  plot(0:nper,mean(ssim(:,:,ip)),'k')
  ytickformat('%.1f')
  title(['Simulated and Expected Wealth: Agent ' num2str(ip)])
  xlabel('Period')
  ylabel('Wealth')
end


%% PARAMETRIC SENSITIVITY ANALYSIS

% Initializations
s1avgplot = zeros(nplot,1);
s1stdplot = zeros(nplot,1);
s2avgplot = zeros(nplot,1);
s2stdplot = zeros(nplot,1);
vinit = v;
xinit = x;

fprintf('Varying Risk Sharing Rate\n')
modelplot = model;
psimin = 0.0;
psimax = 0.2;
psiplot = nodeunif(nplot,psimin,psimax);
v = vinit; x = xinit;
for ip=1:nplot
  modelplot.params = {alpha beta gamma psiplot(ip)};
  [v,x,c,sr,vr,xr] = gamesolve(modelplot,basis,v,x);
  nr  = n*3;
  sr1 = reshape(sr(:,1),nr);
  sr2 = reshape(sr(:,2),nr);
  ss    = cell(2,1);
  ss{1} = sr1(:,1);
  ss{2} = sr2(1,:)';
  sinit = 7.49*ones(nrep,2);
  [ssim,xsim] = gamesimul(modelplot,sinit,nper,ss,xr,esim);
  s1avgplot(ip) = mean(ssim(:,nper,1));
  s1stdplot(ip) = std(ssim(:,nper,1));
  s2avgplot(ip) = mean(ssim(:,nper,2));
  s2stdplot(ip) = std(ssim(:,nper,2));
end
figure
plot(psiplot,[s1avgplot s2avgplot])
ytickformat('%.1f')
xlabel('Risk Sharing Rate')
ylabel('Ergodic Mean Wealth')
legend('Agent 1','Agent 2')
figure
plot(psiplot,[s1stdplot s2stdplot])
ytickformat('%.2f')
xlabel('Risk Sharing Rate')
ylabel('Ergodic Standard Deviation of Wealth')
legend('Agent 1','Agent 2')


%% SAVE FIGURES
printfigures(mfilename)


%% GAMESOLVE FUNCTION FILE
function [out1,out2,out3] = func(flag,i,s,x,e,alpha,beta,gamma,psi)
ns = size(s,1);
if ns>1
  x = squeeze(x);
end
switch flag
  case 'b'
    xl = zeros(ns,1);
    xu = 0.99*s(:,i);
    out1=xl; out2=xu; out3=[];
  case 'f'
    f   = ((s(:,i)-x(:,i)).^(1-alpha(i)))/(1-alpha(i));
    fx  = -(s(:,i)-x(:,i)).^(-alpha(i));
    fxx = -alpha(i)*(s(:,i)-x(:,i)).^(-alpha(i)-1);
    out1=f; out2=fx; out3=fxx;
  case 'g'
    g1   = gamma(1)*x(:,1) + (1-psi)*e(:,1).*x(:,1).^beta(1) + psi*e(:,2).*x(:,2).^beta(2);
    g2   = gamma(2)*x(:,2) + (1-psi)*e(:,2).*x(:,2).^beta(2) + psi*e(:,1).*x(:,1).^beta(1);
    if i==1
      dg1  = gamma(1) + beta(1)*(1-psi)*e(:,1).*x(:,1).^(beta(1)-1);
      dg2  = beta(1)*psi*e(:,1).*x(:,1).^(beta(1)-1);
      ddg1 = (beta(1)-1)*beta(1)*(1-psi)*e(:,1).*x(:,1).^(beta(1)-2);
      ddg2 = (beta(1)-1)*beta(1)*psi*e(:,1).*x(:,1).^(beta(1)-2);
    else
      dg1  = beta(2)*psi*e(:,2).*x(:,2).^(beta(2)-1);
      dg2  = gamma(2) + beta(2)*(1-psi)*e(:,2).*x(:,2).^(beta(2)-1);
      ddg1 = (beta(2)-1)*beta(2)*psi*e(:,2).*x(:,2).^(beta(2)-2);
      ddg2 = (beta(2)-1)*beta(2)*(1-psi)*e(:,2).*x(:,2).^(beta(2)-2);
    end
    g    = [g1   g2];
    gx   = [dg1  dg2];
    gxx  = [ddg1 ddg2];
    gx  = reshape(gx,ns,2,1);
    gxx = reshape(gxx,ns,2,1,1);
    out1=g; out2=gx; out3=gxx;
end
end