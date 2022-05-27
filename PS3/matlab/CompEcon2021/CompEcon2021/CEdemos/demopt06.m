%% DEMOPT06 KKT Conditions for Constrained Optimization Problems

% Preliminary tasks
deminit(mfilename)


x=nodeunif(100,-0.5,1.5);
a=0.1;
b=0.9;

figure

% Plot internal maximum
subplot(1,2,1)
f = @(x) 1.5-2*(x-.75).^2;
hold on
plot(x,f(x))
plot([a;a],[-0.5;2],'k')
plot([b;b],[-0.5;2],'k')
xticks([a b])
set(gca,'XTickLabel',{'a' 'b'})
yticks([])
xstar = 0.75;
ystar = f(xstar);
plotbullet(xstar,ystar)
axis([a b 0.5 2.0],'square')
title('Internal Maximum')

% Plot corner maximum
subplot(1,2,2)
f = @(x) 1-(x+0.25).^2;
hold on
plot(x,f(x))
plot([a;a],[-0.5;2],'k')
plot([b;b],[-0.5;2],'k')
xticks([a b])
set(gca,'XTickLabel',{'a' 'b'})
yticks([])
xstar = a;
ystar = f(xstar);
plotbullet(xstar,ystar)
axis([a b -0.5 1.0],'square')
title('Corner Maximum')


%% SAVE FIGURES
printfigures(mfilename)