%% DEMOPT01 Maximization via Golden Search

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Function to be mazimized and domain
f = @(x) x.*cos(x.^2);
xmin = 0;
xmax = 3;


%% BASE CASE MODEL SOLUTION

% Maximize f from distinct starting values
x1 = golden(f,0,1);
x2 = golden(f,2,3);

% Generate nodes for plotting
x = nodeunif(500,xmin,xmax);

% Plot function and local maxima found by golden search
figure
hold on
plot(x,f(x))
plotbullet(x1,f(x1),20)
plotbullet(x2,f(x2),20)
plottext(x1,f(x1),'Local Max 1','center','bottom',14)
plottext(x2,f(x2),'Local Max 2','center','bottom',14)
xtickformat('%.1f')
title('Maximization of $f(x)=x\cos(x^2)$ via Golden Search')
xlabel('$x$')
ylabel('$f(x)$')


%% SAVE FIGURES
printfigures(mfilename)