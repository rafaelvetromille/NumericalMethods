%% DEMDDP02 Asset Replacement Model

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Base case model parameters
A     = 6;                  % maximum asset age 
alpha = [50 -2.5 -2.5];     % production function coefficients
kappa = 40;                 % net replacement cost
pbar  = 1;                  % long-run mean profit contribution
delta = 0.9;                % discount factor 

% State space
S = (1:A)';                 % asset age
n = length(S);              % total number of states

% Action space (no action=1, replace=2)
X = (1:2)';                 % vector of actions
m = length(X);              % number of actions

% Reward function
f = zeros(n,m);
f(:,1) = pbar*(alpha(1)+alpha(2)*S+alpha(3)*S.^2);
f(:,2) = (pbar*alpha(1)-kappa)*ones(n,1);

% State transition function
g = ones(n,m);
for i=1:n-1
  g(i,1) = i+1;
end

% Model structure
clear model
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% BASE CASE MODEL SOLUTION

% Solve Bellman equation
[v,x,pstar] = ddpsolve(model);
   

%% BASE CASE MODEL SIMULATION

% Simulation parameters
nyrs = 12;                              % number of years simulated

% Initialize simulation
sinit = 1;

% Simulate Model
[spath,xpath] = ddpsimul(pstar,sinit,nyrs,x);

% Plot simulated state path (Age)
figure
plot(0:nyrs,S(spath))
xlim([0 12])
yticks(1:3)
title('Simulated Asset Age')
xlabel('Year')
ylabel('Age of Asset')


%% SAVE FIGURES
printfigures(mfilename)