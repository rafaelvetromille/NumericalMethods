%% Formal beggining of the document
clearvars;
close all;
clc;

%% Exercise 1.

% Calibration
rho      = 0.95;                         
sig      = 0.007;         

N        = 9;        % Number of points (states);  
m        = 3;        % Scaling parameter (I don't know why!); 

% Tauchen's Method

% FUNCTION: [S, P] = tauchen(m, rho, sig, N)
%
% INPUTS:
%     m:    scalar, scaling parameter
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%     n:    scalar, number of grid points for the discretized process
%
% OUTPUTS:
%    S:     column vector of size Nx1, contains all possible states in ascending order
%    P:     matrix of size NxN, contains the transition proabilities. 
%           Rows are current state and columns future state

[S9T, P9T] = tauchen(m, rho, sig, N)

%% Exercise 2.

% Calibration
rho      = 0.95;                         
sig      = 0.007;         

% Number of points (states)
N        = 9;          

% Rouwenhorst's Method

% FUNCTION: [S, Pi] = rouwenhorst(m, rho, sig, N)
% At the end of the document.
%
% INPUTS:
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%   n:      scalar, number of grid points for the discretized process
%
% OUTPUTS:
%   S:      column vector of size Nx1, contains all possible states in ascending order
%   Pi:     matrix of size NxN, contains the transition proabilities. 
%           Rows are current state and columns future state

[S9R, P9R] = rouwenhorst(rho, sig, N)

%% Exercise 3.

%%%% Simulation of the AR(1) process %%%%

% Setting seed, for replicate the results
rng(1234)

% Number of periods 
T = 10000;

% Simulation AR(1) process and shocks

% FUNCTION: [eps, Z] = ar1_simulation(rho, sig, N)
% At the end of the document.
%
% INPUTS:
%   rho:    scalar, AR-coefficient
%   sig:    scalar, standard deviation of innovations
%   N:      scalar, number of periods of time
%
% OUTPUTS:
%   eps:   column vector of size Nx1, contains all shocks
%   Z:     column vector of size Nx1, contains all realizations of Z(t)

[eps, Z] = ar1_simulation(rho, sig, T);

%%%% Simulation of the discretized AR(1) process %%%%

% Initial state
th0 = find(S9T == median(S9T(:)));

% Correspondents indices
IDXT = shock(th0, sig, eps, P9T, T);
IDXR = shock(th0, sig, eps, P9R, T);

% Discretized AR(1) - rho = 0.95
AR1D_T95 = S9T(IDXT);
AR1D_R95 = S9R(IDXR);

% Grid for x
X = 1:T;

% Merge 4 graphics
tiledlayout(2,2)

% Plot Tauchen's Method
nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 1], 'LineWidth', 1.5, 'DisplayName','cos(3x)')
ylim([-0.08 0.08])
%line([0; T].*ones(size(S9T)), [S9T;S9T], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title('Realização do Processo Contínuo para Z(t), \rho = 0.95')
legend('Realização do Processo Contínuo para Z(t)', 'Location', 'southoutside')

nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 0.1], 'LineWidth', 1.5, 'DisplayName','cos(3x)')
hold on
plot(X, AR1D_T95, 'Color', [1, 0.12, 0.12, 1], 'LineWidth', 1, 'DisplayName','cos(3x)')
ylim([-0.08 0.08])
line([0; T].*ones(size(S9T)), [S9T;S9T], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
legend('Realização do Processo Contínuo para Z(t)','Processo Z(t) discretizado via Método de Tauchen', 'Location', 'southoutside')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title({'Realização do Processo Contínuo para Z(t)', 'vs. Discretização via Método de Tauchen'})
hold off 

% Plot Rouwenhorst's Method
nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 1], 'LineWidth', 1.5, 'DisplayName','cos(3x)')
ylim([-0.08 0.08])
%line([0; T].*ones(size(S9R)), [S9R;S9R], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title('Realização do Processo Contínuo para Z(t), \rho = 0.95')
legend('Realização do Processo Contínuo para Z(t)', 'Location', 'southoutside')

nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 0.1], 'LineWidth', 2, 'DisplayName','cos(3x)')
hold on
plot(X, AR1D_R95, 'Color', [1, 0.12, 0.12, 1], 'LineWidth', 1, 'DisplayName','cos(3x)')
ylim([-0.08 0.08])
line([0; T].*ones(size(S9R)), [S9R;S9R], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
legend('Realização do Processo Contínuo para Z(t)','Processo Z(t) discretizado via Método de Rouwenhorst', 'Location', 'southoutside')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title({'Realização do Processo Contínuo para Z(t)', 'vs. Discretização via Método de Rouwenhorst'})
hold off 

%% Exercise 4.

% For simplicity, let's call the variables as Y1 (Tauchen) and Y2 (Rouwenhorst)
Y1 = AR1D_T95;
Y2 = AR1D_R95;

% Run the regressions
lm1_95 = fitlm(lagmatrix(Y1,1), Y1, 'Intercept', false);  % Tauchen, rho = 0.95
lm2_95 = fitlm(lagmatrix(Y2,1), Y2, 'Intercept', false);  % Rouwenhorst, rho = 0.95

%% Exercise 5.1

clearvars;
clc;

% Re-calibration
rho      = 0.99;                         
sig      = 0.007;         

N        = 9;        % Number of points (states);  
m        = 3;        % Scaling parameter (I don't know why!); 

% Tauchen's method
[S9T, P9T] = tauchen(m, rho, sig, N)

%% Exercise 5.2

% Rouwenhorst's method
[S9R, P9R] = rouwenhorst(rho, sig, N)

%% Exercise 5.3

%%%% Simulation of the AR(1) process, rho = 0.99 %%%%

% Setting seed, for replicate the results
rng(1234)

% Number of periods 
T = 10000;

% Simulation AR(1) process and shocks
[Shock, Z] = ar1_simulation(rho, sig, T);

%%%% Simulation of the discretized AR(1) process, rho = 0.99 %%%%

% Initial state (median)
th0 = find(S9T == median(S9T(:)));

% Correspondents indices
IDXT = shock(th0, sig, Shock, P9T, T);
IDXR = shock(th0, sig, Shock, P9R, T);

% Discretized AR(1)
AR1D_T99 = S9T(IDXT);
AR1D_R99 = S9R(IDXR);

% Grid for x
X = 1:T;

% Merge 4 graphics
tiledlayout(2,2)

% Plot Tauchen's Method
nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 1], 'LineWidth', 1.5)
%line([0; T].*ones(size(S9T)), [S9T;S9T], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title('Realização do Processo Contínuo para Z(t), \rho = 0.99')
legend('Realização do Processo Contínuo para Z(t)', 'Location', 'southoutside')

nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 0.1], 'LineWidth', 1.5, 'DisplayName','cos(3x)')
hold on
plot(X, AR1D_T99, 'Color', [1, 0.12, 0.12, 1], 'LineWidth', 2, 'DisplayName','cos(3x)')
line([0; T].*ones(size(S9T)), [S9T;S9T], 'LineWidth', 0.1, 'Color', 'k', 'LineStyle', '--')
legend('Realização do Processo Contínuo para Z(t)','Processo Z(t) discretizado via Método de Tauchen', 'Location', 'southoutside')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title({'Realização do Processo Contínuo para Z(t)', 'vs. Discretização via Método de Tauchen'})
hold off 

% Plot Rouwenhorst's Method
nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 1], 'LineWidth', 1.5, 'DisplayName','cos(3x)')
%line([0; T].*ones(size(S9R)), [S9R;S9R], 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title('Realização do Processo Contínuo para Z(t), \rho = 0.99')
legend('Realização do Processo Contínuo para Z(t)', 'Location', 'southoutside')

nexttile
plot(X, Z, 'Color', [0.128, 0.128, 0.128, 0.1], 'LineWidth', 1, 'DisplayName','cos(3x)')
hold on
plot(X, AR1D_R99, 'Color', [1, 0.12, 0.12, 1], 'LineWidth', 2, 'DisplayName','cos(3x)')
line([0; T].*ones(size(S9R)), [S9R;S9R], 'LineWidth', 0.1, 'Color', 'k', 'LineStyle', '--')
legend('Realização do Processo Contínuo para Z(t)','Processo Z(t) discretizado via Método de Rouwenhorst', 'Location', 'southoutside')
xlabel('Período de Tempo')
ylabel('Realização de Z(t)')
title({'Realização do Processo Contínuo para Z(t)', 'vs. Discretização via Método de Rouwenhorst'})
hold off 

%% Exercise 5.4

% For simplicity, let's call the variables as Y1 (Tauchen) and Y2 (Rouwenhorst)
Y1 = AR1D_T99;
Y2 = AR1D_R99;

% Run the regressions
lm1_99 = fitlm(lagmatrix(Y1,1), Y1, 'Intercept', false);
lm2_99 = fitlm(lagmatrix(Y2,1), Y2, 'Intercept', false);

%% FUNCTIONS 

% tauchen.m
function [theta, Pi] = tauchen(m, rho, sig, N)

% Setting the grid and the trasition matrix
theta   = zeros(N,1)';
Pi      = eye(N,N);

% Specifying the grid
theta(N)    = + m * sqrt(sig^2/(1-rho^2)); 
theta(1)    = - m * sqrt(sig^2/(1-rho^2)); 

% Computing the remaining (N-2) points [equidistantly distributed]  
step     = (theta(N)-theta(1))/(N-1);

for i = 2:(N-1)
   theta(i) = theta(i-1) + step;
end

% Computing the transition probabilities
for j = 1:N
    for k = 1:N
        % The transition to the corner points have to be treated differently
        if k == 1
            Pi(j,k) = normcdf((theta(1) - rho*theta(j) + step/2) / sig);
        % The transition to the corner points have to be treated differently
        elseif k == N
            Pi(j,k) = 1 - normcdf((theta(N) - rho*theta(j) - step/2) / sig);
        else
        % The other points will be calculated by:
           Pi(j,k) = normcdf((theta(k) - rho*theta(j) + step/2) / sig) - ...
			   normcdf((theta(k) - rho*theta(j) - step/2) / sig);
        end
    end
end
end

% Rouwenhorst.m
function [theta, Pi] = rouwenhorst(rho, sig, N)

% Specifying the grid
theta(N)    = + (sig / sqrt(1-rho^2)) * sqrt(N-1); 
theta(1)    = - (sig / sqrt(1-rho^2)) * sqrt(N-1);

% Computing the remaining (N-2) points [equidistantly distributed]  
step        = (theta(N)-theta(1))/(N-1);

for i = 2:(N-1)
   theta(i) = theta(i-1) + step;
end

% Computing p and P2
p  = (1+rho)/2;
Pi = [p 1-p; 1-p p];

% Compute the transition matrix P recursively
for n = 3:N
    Pi =    p * [Pi zeros(n-1,1); zeros(1,n-1) zeros(1,1)] + ...
        (1-p) * [zeros(n-1,1) Pi; zeros(1,1) zeros(1,n-1)] + ...
        (1-p) * [zeros(1,n-1) zeros(1,1); Pi zeros(n-1,1)] + ...
            p * [zeros(1,1) zeros(1,n-1); zeros(n-1,1) Pi]; 
end

% Normalize rows to add up to 1
for i = 1:N, Pi(i,:) = Pi(i,:)/sum(Pi(i,:)); end
end

% ar1_simulation.m
function [eps, Z] = ar1_simulation(rho, sig, N)
    
    Z   = zeros(N,1)';
    eps = normrnd(0, sig, N, 1)';
    
    for i=2:N
        Z(i) = rho*Z(i-1) + eps(i);
    end
end

% shock.m
function idx = shock(th0, sig, eps, Pi, T)

idx    = ones(T,1);
idx(1) = th0;
cum    = cumsum(Pi')';

for i = 2:T
    x     = find(normcdf(eps(i), 0, sig) <= cum(idx(i-1),:));
    idx(i) = x(1);
end
end