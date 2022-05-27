%% DEMMATH03 Discrete and Continuous Distributions
%
% Draws probability density, probability mass, and cumulative distribution
% functions for common discrete and continuous random variables.
% Preliminary tasks

% Preliminary tasks
deminit(mfilename)


%% Standard Normal Distribution

% Distribution parameters
mu  = 0;
var = 1;

% Plotting grid
n = 400;
xmin = -4;
xmax =  4;
x = nodeunif(n,xmin,xmax);

% Compute pdf, cdf, and 90% quantile
f = pdf('Normal',x,mu,var);
F = cdf('Normal',x,mu,var);
xquant = icdf('Normal',0.9,mu,var);

% Plot pdf
figure
hold on
plot(x,f)
ytickformat('%.1f')
title('Standard Normal Probability Density Function')
xlabel('$x$')
ylabel('Density')

% Plot cdf
figure
hold on
plot(x,F)
plotvdash(xquant,0.9)
plothdash(xquant,0.9)
plottext(-4,0.9,'$0.9$','left','bottom',14)
plottext(xquant,0,'$x_{0.9}$','left','bottom',16)
ytickformat('%.1f')
ylabel('Probability')
title('Standard Normal Cumulative Distribution Function')
xlabel('$x$')


%% Lognormal Distribution

% Distribution parameters
mu1  = 0.0;
var1 = 0.5;
mu2  = 0.0;
var2 = 1.0;
mu3  = 0.5;
var3 = 1.0;

% Plotting grid
n = 400;
xmin = 0;
xmax = 6;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Lognormal',x,mu1,var1);
F1 = cdf('Lognormal',x,mu1,var1);
f2 = pdf('Lognormal',x,mu2,var2);
F2 = cdf('Lognormal',x,mu2,var2);
f3 = pdf('Lognormal',x,mu3,var3);
F3 = cdf('Lognormal',x,mu3,var3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
ytickformat('%.1f')
title('Lognormal Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$\mu=0.0,\sigma^2=0.5$','$\mu=0.0,\sigma^2=1.0$','$\mu=0.5,\sigma^2=1.0$')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
ytickformat('%.1f')
xlabel('$x$')
ylabel('Probability')
title('Lognormal Cumulative Distribution Function')
legend('$\mu=0.0,\sigma^2=0.5$','$\mu=0.0,\sigma^2=1.0$','$\mu=0.5,\sigma^2=1.0$')


%% Exponential Distribution

% Distribution parameters
beta1 = 1;
beta2 = 2;

% Plotting grid
n = 400;
xmin = 0;
xmax = 5;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Exponential',x,beta1);
F1 = cdf('Exponential',x,beta1);
f2 = pdf('Exponential',x,beta2);
F2 = cdf('Exponential',x,beta2);

% Plot pdf
figure
hold on
plot(x,[f1 f2])
ytickformat('%.1f')
title('Exponential Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$\beta=1$','$\beta=2$')

% Plot cdf
figure
hold on
plot(x,[F1 F2])
ytickformat('%.1f')
title('Exponential Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$\beta=1$','$\beta=2$')

 
%% Gamma Distribution

% Distribution parameters
alpha1 = 1.0;
theta1 = 0.5;
alpha2 = 2.0;
theta2 = 0.5;

% Plotting grid
n = 400;
xmin = 0;
xmax = 4;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Gamma',x,alpha1,theta1);
F1 = cdf('Gamma',x,alpha1,theta1);
f2 = pdf('Gamma',x,alpha2,theta2);
F2 = cdf('Gamma',x,alpha2,theta2);

% Plot pdf
figure
hold on
plot(x,[f1 f2])
ytickformat('%.1f')
title('Gamma Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$\alpha=1,\theta=0.5$','$\alpha=2,\theta=0.5$')

% Plot cdf
figure
hold on
plot(x,[F1 F2])
ytickformat('%.1f')
title('Gamma Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$\alpha=1,\theta=0.5$','$\alpha=2,\theta=0.5$')


%% Beta Distribution

% Distribution parameters
alpha1 = 3;
beta1  = 3;
alpha2 = 2;
beta2  = 5;
alpha3 = 5;
beta3  = 2;

% Plotting grid
n = 400;
xmin = 0;
xmax = 1;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Beta',x,alpha1,beta1);
F1 = cdf('Beta',x,alpha1,beta1);
f2 = pdf('Beta',x,alpha2,beta2);
F2 = cdf('Beta',x,alpha2,beta2);
f3 = pdf('Beta',x,alpha3,beta3);
F3 = cdf('Beta',x,alpha3,beta3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
xtickformat('%.1f')
ytickformat('%.1f')
ylabel('Density')
title('Beta Probability Density Function')
xlabel('$x$')
legend('$\alpha=3,\beta=3$','$\alpha=2,\beta=5$','$\alpha=5,\beta=2$')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
xtickformat('%.1f')
ytickformat('%.1f')
title('Beta Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$\alpha=3,\beta=3$','$\alpha=2,\beta=5$','$\alpha=5,\beta=2$')


%% Chi-Squared Distribution

% Distribution parameters
k1 = 2;
k2 = 3;
k3 = 5;

% Plotting grid
n = 400;
xmin = 0;
xmax = 12;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Chisquare',x,k1);
F1 = cdf('Chisquare',x,k1);
f2 = pdf('Chisquare',x,k2);
F2 = cdf('Chisquare',x,k2);
f3 = pdf('Chisquare',x,k3);
F3 = cdf('Chisquare',x,k3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
ytickformat('%.1f')
title('Chi-Squared Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$k=2$','$k=3$','$k=5$')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
ytickformat('%.1f')
title('Chi-Squared Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$k=2$','$k=3$','$k=5$')


%% F Distribution

% Distribution parameters
n1 = 2;
m1 = 1;
n2 = 10;
m2 = 1;
n3 = 100;
m3 = 100;

% Plotting grid
n = 400;
xmin = 0;
xmax = 5;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('F',x,n1,m1);
F1 = cdf('F',x,n1,m1);
f2 = pdf('F',x,n2,m2);
F2 = cdf('F',x,n2,m2);
f3 = pdf('F',x,n3,m3);
F3 = cdf('F',x,n3,m3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
ytickformat('%.1f')
title('F Probability Density Function')
legend('$n=2,m=1$','$n=10,m=1$','$n=100,m=100$')
xlabel('$x$')
ylabel('Density')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
ytickformat('%.1f')
title('F Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$n=2,m=1$','$n=10,m=1$','$n=100,m=100$')


%% Student's T Distribution

% Distribution parameters
k1 = 1;
k2 = 4;
k3 = 100;

% Plotting grid
n = 400;
xmin = -4;
xmax =  4;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('T',x,k1);
F1 = cdf('T',x,k1);
f2 = pdf('T',x,k2);
F2 = cdf('T',x,k2);
f3 = pdf('T',x,k3);
F3 = cdf('T',x,k3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
ytickformat('%.1f')
title('Student''s T Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$\nu=1$','$\nu=4$','$\nu=100$')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
ytickformat('%.1f')
title('Student''s T Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$\nu=1$','$\nu=4$','$\nu=100$')


%% Logistic Distribution

% Distribution parameters
mu1    = 0;
sigma1 = 1;
mu2    = 2;
sigma2 = 1;
mu3    = 2;
sigma3 = 2;

% Plotting grid
n = 400;
xmin = -15;
xmax =  15;
x = nodeunif(n,xmin,xmax);

% Compute pdf & cdf
f1 = pdf('Logistic',x,mu1,sigma1);
F1 = cdf('Logistic',x,mu1,sigma1);
f2 = pdf('Logistic',x,mu2,sigma2);
F2 = cdf('Logistic',x,mu2,sigma2);
f3 = pdf('Logistic',x,mu3,sigma3);
F3 = cdf('Logistic',x,mu3,sigma3);

% Plot pdf
figure
hold on
plot(x,[f1 f2 f3])
ytickformat('%.2f')
title('Logistic Probability Density Function')
xlabel('$x$')
ylabel('Density')
legend('$\mu=0,\sigma=1$','$\mu=2,\sigma=1$','$\mu=2,\sigma=2$')

% Plot cdf
figure
hold on
plot(x,[F1 F2 F3])
ytickformat('%.1f')
title('Logistic Cumulative Distribution Function')
xlabel('$x$')
ylabel('Probability')
legend('$\mu=3,\sigma=3$','$\mu=2,\sigma=5$','$\mu=5,\sigma=2$')


%% Binomial Distribution

% Distribution parameters
n = 3;
p = 0.6;

% Plotting grid
x = 0:n; 

% Probability mass function
f = pdf('Binomial',x,n,p);

% Plot pmf
figure
hold on
bar(x,f,0.8)
set(gca,'XTick',0:n)
ylim([0 0.5])
ytickformat('%.1f')
title('Binomial Probability Mass Function')
xlabel('$x$')
ylabel('Mass')


%% Geometric Distribution

% Distribution parameters
n = 8;
p = 0.4;

% Plotting grid
x = 0:n; 

% Probability mass function
f = pdf('Geometric',x,p);

% Plot pmf
figure
hold on
bar(x,f,0.8)
ytickformat('%.2f')
title('Geometric Probability Mass Function')
xlabel('$x$')
ylabel('Mass')


%% Poisson Distribution

% Distribution parameters
n = 7;
lambda = 1.5;

% Plotting grid
x = 0:n; 

% Probability mass function
f = pdf('Poisson',x,lambda);

% Plot pmf
figure
hold on
bar(x,f,0.8)
set(gca,'XTick',0:n)
ytickformat('%.2f')
title('Poisson Probability Mass Function')
xlabel('$x$')
ylabel('Mass')


%% SAVE FIGURES
printfigures(mfilename)