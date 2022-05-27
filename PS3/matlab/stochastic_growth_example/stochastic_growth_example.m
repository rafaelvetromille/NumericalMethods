%%%%% solving the stochastic growth model with collocation methods
clearvars;
close all;
clc;

%%

%%%%% economic parameters

alpha = 1/3;        %% capital's share in production function
beta  = 0.95;       %% time discount factor   
delta = 0.05;       %% depreciation rate
sigma = 1;          %% CRRA (=1/IES)

rhoz  = 0.95;       %% AR1 coefficient, productivity shocks
sigz  = 0.1;        %% innovation std dev, productivity shocks

%%%%% numerical parameters

max_iter = 500;      %% maximum number of iterations
tol      = 1e-7;     %% treat numbers smaller than this as zero

%%%%% grid for productivity shocks

nz       = 7;       %% number of breakpoints for z grid
nez      = 3;       %% number of nodes for quadrature (shocks)

% quadrature nodes and weights
[ez, wz] = qnwunif(nez,1e-9,1-1e-9);
ez       = sigz*norminv(ez,0,1); 

zmin     = exp(-4*sqrt(1/(1-rhoz^2))*sigz);
zmax     = exp( 4*sqrt(1/(1-rhoz^2))*sigz);
zgrid    = exp(nodeunif(nz, log(zmin), log(zmax)));

%%

%%%%% grid for capital stock

nk       = 99;       %% number of breakpoints for k grid
curv     = 0.50;     %% (curv = 0 log-spaced, curv = 1 linear)

kmin     = 1e-3;
kmax     = (zmax/delta)^(1/(1-alpha));
kgrid    = nodeunif(nk, kmin.^curv, kmax.^curv).^(1/curv);

%%%%% put in a structure to pass to other functions

parameters.alpha = alpha;
parameters.beta  = beta;
parameters.delta = delta;
parameters.sigma = sigma;
parameters.rhoz  = rhoz;
parameters.sigz  = sigz;

parameters.ez    = ez;
parameters.wz    = wz;

%%%%% setup state space using CompEcon tools

fspace   = fundef({'spli', kgrid},...
                  {'spli', zgrid}); % function space structure

grid     = funnode(fspace);% nodes where we solve the problem
Phi      = funbas(fspace); % matrix of collocation basis vectors 
                           % Phi_{ij} = phi_j(k_i) 

kgrid    = grid{1}; % extra 2 points for 3rd-order spline
zgrid    = grid{2}; % extra 2 points for 3rd-order spline

kmin     = kgrid(1);
kmax     = kgrid(end);

zmin     = zgrid(1); 
zmax     = zgrid(end); 

%%%%% bounds on consumption

cmin     = tol;
cmax     = zmax*(kmax.^alpha) + (1-delta)*kmax - eps;

%%%%% form collection of states

s        = gridmake(grid); % ns-by-2 matrix where ns=nk*nz
ns       = size(s,1);        

k        = s(:,1);
z        = s(:,2);

%%%%% initial guess at collocation coefficients "a"

c = alpha*beta*z.*k.^alpha; % guess for consumption policy
v = log(c)/(1-beta);        % guess for value function

a = Phi\v;                  % implied collocation coefficients

%%

tic
%%%%% solve Bellman equation

for i=1:max_iter

%%%%% optimal consumption given these coefficints

c = solve_brent('rhs_bellman',s,parameters,a,fspace,cmin,cmax,tol);

%%%%% maximized rhs of Bellman equation

v = rhs_bellman(c,s,parameters,a,fspace); %% v(a)

%%%%% update collocation coefficients
Jacobian = 0;

%%%%% implied by optimal consumption
kprime = z.*(k.^alpha) + (1-delta)*k-c;

%%%%% build up Jacobian matrix of v(a)
for j=1:numel(ez)

zprime = max(min(z.^rhoz.*exp(ez(j)), zmax), zmin);

sprime = [kprime,zprime];

Jacobian = Jacobian+beta*wz(j)*funbas(fspace,sprime);

end

%%%%% Newton's method
anew = a - (Phi-Jacobian)\(Phi*a-v); 

%%%%% check if converged

error = norm(anew-a,inf);

fprintf('%4i %6.2e \n',[i, error]);

if error<tol, break, end;

%%%%% if not converged, update and try again

a = anew;

end

time      = toc;
last_iter = i;

if error<tol && last_iter<=max_iter
   
fprintf('\n');    
fprintf('time taken (seconds) = %3.3f \n',time);
fprintf('number of iterations = %3.0f \n',last_iter);
fprintf('error at last step   = %3.10f \n',error);
fprintf('\n'); 

end

%%

%%%%% reshape solution

Nk = numel(kgrid);
Nz = numel(zgrid);

VV = reshape(v,Nk,Nz);      %% VV = v(k_i,z_j)
CC = reshape(c,Nk,Nz);      %% CC = c(k_i,z_j)

%%%%% optimal policy k'=g(k,z)
GG = reshape(kprime,Nk,Nz); %% GG = g(k_i,z_j)

%%

figure(1)
plot(kgrid,VV,'b-')
xlabel('capital stock, k')
ylabel('value function, v(k,z)')
axis([0 600 -30 60])

figure(2)
plot(kgrid,GG,'r-',k,k,'k--')
xlabel('capital stock, k')
ylabel('policy function, g(k,z)')
axis([0 60 0 60])
