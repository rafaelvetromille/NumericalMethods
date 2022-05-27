%% DEMOPT03 Nelder-Mead Simplex Method
%
% Creates and plays a movie of Nelder-Meade simplex iterations when
% maximizing banana function f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1].

% Preliminary tasks
deminit(mfilename)
warning off


%  The "banana" or Rosencrantz function
banana = @(x) -100*(x(2,:)-x(1,:).^2).^2-(1-x(1,:)).^2;

n = [20 20];
xmin = [-0.2 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);

y = banana(x');
y = reshape(y,n)';
conts = -exp(0.25:0.5:20);

figure
hold on
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlim([-0.21 1.21])
ylim([-0.21 1.21])
plotbullet(1,1,14,'k')
xticks(-0.2:0.2:1.2)
xtickformat('%.1f')
ytickformat('%.1f')
title('Nelder-Mead Maximizes the Banana Function')
xlabel('$x_1$')
ylabel('$x_2$')

optset('neldmead','maxit',1)
x = [1;0];
[xx,S] = neldmead(banana,x);
hp = patch(S(1,:),S(2,:),[0.5 0.5 0.5]);
for i=1:60
  xvec(:,i) = x;
  [x,S] = neldmead(banana,x,S);
  set(hp,'xdata',S(1,:)','ydata',S(2,:)');
  getframe;
  pause(0.1)
end