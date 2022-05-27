%% ODEPHASE
%
%  Plots two-dimensional phase diagram for ODE
%
%  Usage
%    odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,xspx,xstst,x1null,x2null,nomovie)
%  Let
%    d  = dimension of state process x
%    k  = number of initial states
%  Input
%    x1lim    : lower & upper bound of x1 axis
%    x2lim    : lower & upper bound of x2 axis
%    figtitle : figure title
%    x1label  : x1-axis label
%    x2label  : x2-axis label
%    f        : velocity function
%    x        : N.d.k solution values
%    xspx     : N.d.1 separatrix
%    xstst    : d.k steady-state
%    x1null   : N.d x1 nulcline
%    x2null   : N.d x2 nulcline
%    nomovie  : 1 if no movie, empty otherwise
%  Output
%    two-dimensional phase diagram
%
%  Velocity Function
%    User-supplied function that returns velocity at states according to
%    Format
%      v = f(x)
%    Input
%      x      : d.k states
%    Output
%      v      : d.k velocities
%
%  Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu

function odephase(x1lim,x2lim,ftitle,x1label,x2label,f,x,xspx,xstst,x1null,x2null,nomovie)

figure
hold on
xlim(x1lim)
ylim(x2lim)
if ~isempty(ftitle),  title(ftitle),   end
if ~isempty(x1label), xlabel(x1label), end
if ~isempty(x2label), ylabel(x2label), end

if nargin<6
  return
end

% Plot Velocity Field
xg = gridmake(nodeunif(11,x1lim(1),x1lim(2)),nodeunif(11,x2lim(1),x2lim(2)))';
v = real(f(xg));
quiver(xg(1,:),xg(2,:),v(1,:),v(2,:),'o','LineWidth',1,'MarkerSize',3,'Color',[0.6 0.6 0.6]);

if nargin<7
  return
end

% Plot Separatrix
if ~isempty(xspx)
  ps = plot(xspx(:,1),xspx(:,2),'k--','LineWidth',2);
end

% Plot Nullclines
p1 = plot(x1null(:,1),x1null(:,2),'c','LineWidth',2);
p2 = plot(x2null(:,1),x2null(:,2),'g','LineWidth',2);

% Plot State Path
if ~isempty(x)
  n = size(x,1);
  m = max(1,floor(n/100));
  for k=1:size(x,3)
    plotbullet(x(1,1,k),x(1,2,k),18,'r')
    if nargin<12
      for j=m+1:n
        if mod(j-1,m)==0
          plot(x(j-m:j,1,k),x(j-m:j,2,k),'r','LineWidth',2)
          getframe;
          if x(j-m,1,k)<x1lim(1)||x(j-m,1,k)>x1lim(2)||x(j-m,2,k)<x2lim(1)||x(j-m,2,k)>x2lim(2), break, end
        end
      end
    else
      plot(x(:,1,k),x(:,2,k),'r','LineWidth',2)
    end
  end
end

% Plot Steady State
if ~isempty(xstst)
  for k=1:size(xstst,2)
    plotbullet(xstst(1,k),xstst(2,k))
  end
end

% Legend
if ~isempty(xspx)
  legend([ps p1 p2],'Separatrix',[x1label ' ' 'Nullcline'],[x2label ' ' 'Nullcline'])
else
  legend([p1 p2],[x1label ' ' 'Nullcline'],[x2label ' ' 'Nullcline'])
end