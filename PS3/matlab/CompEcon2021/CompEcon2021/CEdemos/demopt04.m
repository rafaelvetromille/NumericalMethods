%% DEMOPT04 Maximization of Rosencrantz Function by Various Methods

% f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1].

% Preliminary tasks
deminit(mfilename)
warning off
optset('qnewton','maxit',1)
optset('qnewton','ShowIters',0)
optset('neldmead','maxit',1)
optset('neldmead','ShowIters',0)

%  Rosencrantz function
f  = @(x) -100*(x(2,:)-x(1,:).^2).^2-(1-x(1,:)).^2;

n = [40 40];
xmin = [-0.7 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);
y = f(x');
y = reshape(y,n(1),n(2))';
conts = -exp(0.25:0.5:20);


%% Steepest Ascent Maximization

optset('qnewton','SearchMeth',1)
k = 250;
x = [1;0];
[~,A] = qnewton(f,x);
for i=1:k
  xx1(:,i) = x;
  [x,A] = qnewton(f,x,A);
  if norm(x-[1;1])>sqrt(eps), iters1=i; end
end

figure
hold on
plot(xx1(1,:),xx1(2,:))
plotbullet(xx1(1,:),xx1(2,:),18,'r')
plotbullet(1,1,18,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xtickformat('%.1f')
ytickformat('%.1f')
yticks([0 0.5 1])
axis([-0.7 1.2 -0.2 1.2])
title('Steepest Ascent Maximization of Banana Function')
xlabel('$x_1$','VerticalAlignment','cap')
ylabel('$x_2$','VerticalAlignment','bottom')


%% DFP Maximization

optset('qnewton','SearchMeth',2)
x = [1;0];
[~,A] = qnewton(f,x);
for i=1:k
  xx2(:,i) = x;
  [x,A] = qnewton(f,x,A);
  if norm(x-[1;1])>sqrt(eps), iters2=i; end
end

figure
hold on
plot(xx2(1,:),xx2(2,:))
plotbullet(xx2(1,:),xx2(2,:),18,'r')
plotbullet(1,1,18,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xtickformat('%.1f')
ytickformat('%.1f')
yticks([0 0.5 1])
axis([-0.7 1.2 -0.2 1.2])
title('DFP Quasi-Newton Maximization of Banana Function')
xlabel('$x_1$')
ylabel('$x_2$')


%% BFGS Maximization

optset('qnewton','SearchMeth',3)
x = [1;0];
[~,A] = qnewton(f,x);
for i=1:k
  xx3(:,i) = x;
  [x,A] = qnewton(f,x,A);
  if norm(x-[1;1])>sqrt(eps), iters3=i; end
end

figure
hold on
plot(xx3(1,:),xx3(2,:))
plotbullet(xx3(1,:),xx3(2,:),18,'r')
plotbullet(1,1,18,'k')
yticks([0 0.5 1])
axis([-0.7 1.2 -0.2 1.2])
contour(xcoord{1},xcoord{2},y,conts,'k:')
xtickformat('%.1f')
ytickformat('%.1f')
title('BFGS Quasi-Newton Maximization of Banana Function')
xlabel('$x_1$')
ylabel('$x_2$')


%% Nelder-Mead Maximization

x = [1;0];
[~,S] = neldmead(f,x);
for i=1:60
  xx4(:,i) = x;
  [x,S] = neldmead(f,x,S);
  if norm(x-[1;1])>sqrt(eps), iters4=i; end
end

figure
hold on
plot(xx4(1,:),xx4(2,:))
plotbullet(xx4(1,:),xx4(2,:),18,'r')
plotbullet(1,1,18,'k')
contour(xcoord{1},xcoord{2},y,conts,'k:')
xtickformat('%.1f')
ytickformat('%.1f')
yticks([0 0.5 1])
axis([-0.7 1.2 -0.2 1.2])
title('Nelder-Mead Maximization of Banana Function')
xlabel('$x_1$')
ylabel('$x_2$')

fprintf('\n')
fprintf('Iterations Required by Method\n');
fprintf('Steepest    %8i \n',iters1)
fprintf('DFP         %8i \n',iters2)
fprintf('BFGS        %8i \n',iters3)
fprintf('Nelder-Mead %8i \n',iters4)


%% SAVE FIGURES
printfigures(mfilename)