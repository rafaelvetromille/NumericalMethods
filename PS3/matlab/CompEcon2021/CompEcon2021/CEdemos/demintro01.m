%% DEMINTRO01 Inverse Demand Problem
%
% Plots demand function q(p)=0.5*p^-0.2+0.5*p^-0.5 and its inverse, and
% computes price that will clear the market of a quantity 2,

% Preliminary tasks
deminit(mfilename)


% Target quantity and plotting limits
qstar = 2;
pmin  = 0.02;
pmax  = 0.40;

% Compute price that clears market of q=2 using Newton's method
p = 0.25;
fprintf('\n')
for it=1:100
  f = 0.5*p^-0.2 + 0.5*p^-0.5 - qstar;
  d = -0.01*p^-1.2 - 0.25*p^-1.5;
  s = -f/d;
  p = p + s;
  fprintf('iteration %3i  price %8.4f\n',[it p])
  if norm(s)<1.e-8, break, end
end
pstar = p;

% Generate demand function
n = 100;
p = nodeunif(n,pmin,pmax);
q = 0.5*p.^-0.2 + 0.5*p.^-0.5;

% Graph demand and inverse demand functions
figure
subplot(1,2,1)
plot(p,q)
xticks([])
yticks([])
axis('square')
box off
title('Demand')
xlabel('$p$')
ylabel('$q$')
subplot(1,2,2)
hold on
plot(q,p)
plotvdash(qstar,pstar)
plothdash(qstar,pstar)
plotbullet(qstar,pstar)
plottext(1.3,pstar,'$p^*$','right','middle',14);
plottext(qstar,0,'$2$','center','top',14);
xticks([])
yticks([])
axis('square')
box off
title('Inverse Demand')
xlabel('$q$')
ylabel('$p$')


%% Save Figures
printfigures(mfilename)