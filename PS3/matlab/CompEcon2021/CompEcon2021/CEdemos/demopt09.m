%% DEMOPT09 Linear Programming lpsolve

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Base case model parameters
A = [3 1; 1 2; 1 1];
b = [90; 80; 40];
p = [30; 20];


%% BASE CASE MODEL SOLUTION

[x,fval,lambda,MP,surplus,exitflag] = lpsolve(p,A,b);

fprintf('\n')
fprintf('Maximum Objective  %7.0f\n',fval)

fprintf('\n')
fprintf('            Marginal\n')
fprintf('      X      Profit\n')
fprintf('%7.0f  %7.0f\n',[x MP]')

fprintf('\n')
fprintf('    Shadow\n')
fprintf('    Price   Surplus\n')
fprintf('%7.0f   %7.0f\n',[lambda surplus]')