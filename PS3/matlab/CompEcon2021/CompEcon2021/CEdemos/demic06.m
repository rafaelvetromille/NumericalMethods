%% DEMIC06 Optimal Fish Harvest Model

% Preliminary tasks
deminit(mfilename)


%% Formulation

% Base case model parameters
alpha = 0.5;
sigma = 0.5;
P     = 1;
c     = 0.25;
rho   = 0.1;

% Model structure
clear model
model.func   = @func;
model.params = {alpha sigma P c rho};
model.xindex = [0 0;2 1];
model.F      = [0;0];

% Define starting values
x = [0.01 0;1 0];


%% BASE CASE MODEL SOLUTION

% Call solver
n = 50;
[cv,basis,x] = icsolve(model,x,n,'cheb');
sstar = x(2,1);

% Plot results
s = linspace(min(0.01,x(1,1)),1.5,101)';
s = sort([s;sstar;sstar+sqrt(eps)]);

V = funeval(cv,basis,s);
Vs = funeval(cv,basis,s,1);
Vss = funeval(cv,basis,s,2);

ind = s>sstar;
V(ind) = P*(s(ind)-sstar)-c*log(s(ind)/sstar)+funeval(cv,basis,sstar);
Vs(ind) = P-c./s(ind);
Vss(ind) = c./s(ind).^2;

figure
plot(s,V,'k',sstar,funeval(cv,basis,sstar),'k*')
xtickformat('%.1f')
ytickformat('%.1f')
title('Value Function')
xlabel('$S$')
ylabel('$V$')

figure
plot(s,P-c./s,'r:',s,Vs,'k',sstar,funeval(cv,basis,sstar,1),'k*')
ylim([0 5])
xtickformat('%.1f')
title('Shadow Price Function')
xlabel('$S$')
ylabel('$V''$')

figure
plot(s,Vss,'k',sstar,funeval(cv,basis,sstar,2),'k*')
ylim([-10 5])
xtickformat('%.1f')
title('Curvature of Value Function')
xlabel('$S$')
ylabel('$V''''$');


%% BASE CASE MODEL SIMULATION

% Long run density of fish stocks
[cp,Ex] = itodensity(model,basis);
p = funeval(cp,basis,s);
p(s>sstar) = 0;
basis0 = fundefn('cheb',51,basis.a,3);
s0 = linspace(basis0.a,basis0.b,201)';
cp0 = itodensity(model,basis0);
p0 = funeval(cp0,basis0,s0);

figure
plot(s,p,s0,p0)
yy = get(gca,'ylim');yy(1) = 0; set(gca,'ylim',yy);
xtickformat('%.1f')
ytickformat('%.1f')
title('Long-run Density of Fish Stocks');
xlabel('$S$')
ylabel('Probability')
legend({'With harvesting','Without harvesting'})

disp('      S*       E[S]')
disp([sstar Ex])


%% SAVE FIGURES
printfigures(mfilename)


%% ICSOLVE FUNCTION FILE
function out = func(flag,s,alpha,sigma,P,c,rho)
switch flag
  case 'f'        % reward
    out = zeros(size(s,1),1);
  case 'mu'       % state drift
    out = alpha*(1-s).*s;
  case 'sigma'    % state diffusion
    if sigma==0
      out = [];
    else
      out = sigma*s;
    end
  case 'rho'      % state contingent discount rate
    out = rho+zeros(size(s,1),1);
  case 'R+'       % reward associated with trigger s(1) and target s(2) (s(1)<s(2))
    out = zeros(1,4);
  case 'R-'       % reward associated with trigger s(1) and target s(2) (s(1)>s(2))
    S0 = s(:,1); S1 = s(:,2);
    out = [P*(S0-S1)-c*log(S0./S1) P-c./S0 -P+c./S1 c./S0.^2];
end
end