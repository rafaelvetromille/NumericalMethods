%% DEMOPT08 Constrained Optimization Using nlpsolve

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

f = @(x) -x(1)^2 - (x(2)-1)^2 - 3*x(1) + 2;
g = @(x) [4*x(1)+x(2);x(1)^2+x(2)*x(1)];
b = [0.5;2.0];
x = [0;1];


%% BASE CASE MODEL SOLUTION

[x,fval,lambda,MP,surplus,exitflag] = nlpsolve(x,f,g,b);

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