%% DEMMATH01 Taylor Approximations
%
% Illustrate uni- and bi-variate Taylor approximations.

% Preliminary tasks
deminit(mfilename)


%% Univariate Taylor Approximation

% Compute first- and second-order approximations
x = nodeunif(100,-1,1);
y = (x+1).*exp(2*x);    % function
y1 = 1+3*x;             % first derivative
y2 = 1+3*x+4*x.^2;      % second derivative

% Plot first- and second-order approximations
figure
plot(x,[y y1 y2])
xtickformat('%.1f')
title('Taylor Approximations for Univariate Fuction')
xlabel('$x$')
ylabel('$y$')
legend('Function','1st Order Approximation','2nd Order Approximation')


%% Bivariate Taylor Approximation

% Plot function
[x1,x2] = ndgrid(0:0.1:2,-1:0.01:1);
f  = (x1.^2).*exp(-x2);
figure
mesh(x1,x2,f) 
title('$f(x_1,x_2) = x_1^2 - 2x_1x_2 + 0.5x_2^2 + x_2$')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$y$') 

% Compute first- and second-order approximations
[x1,x2] = ndgrid(0.99:0.0001:1.01,-0.01:0.0001:0.01);
f  = (x1.^2).*exp(-x2);
f1 = 2*x1 - x2 - 1;
f2 = x1.^2 - 2*x1.*x2 + 0.5*x2.^2 + x2;

% Plot first-order approximation error
figure
mesh(x1,x2,f1-f) 
xtickformat('%.2f')
ytickformat('%.2f')
title('First-Order Approximation Error')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$y$') 

% Plot second-order approximation error
figure
mesh(x1,x2,f2-f)
xtickformat('%.2f')
ytickformat('%.2f')
title('Second-Order Approximation Error')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$y$')


%% SAVE FIGURES
printfigures(mfilename)