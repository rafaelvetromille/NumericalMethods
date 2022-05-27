%% ODESPX
%
%  Generates separatrix through saddle point of 2-dimensional ODE
%
%  Specifically, generates velocity field for 2-dimensionnal 1st-order ODE
%    x'(t) = f(x(t)), t in [0,T]
%  by solving the ODE backwards in time, starting from near the saddle
%  point in the directions of the stable eigenvector. Here, x is 2.1
%  vector-valued function defined on time domain [0,T] and x' is its 2.1
%  vector-valued derivative with respect to time.
%
%  Usage
%    x = odespx(f,x,T)
%  Input
%    f        : velocity function (see below)
%    x        : 2.1 saddle point
%    T        : time horizon in direction of stable eigenvector
%  Output
%    x        : n.2 separatrix
%  NAN and infinite values removed from x before returning.

%  Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function x = odespx(f,x0,T)

n = 1000;

J = fdjac(f,x0);
[V,D] = eig(J);
i = find(diag(D)<0);
j = find(diag(D)>0);
if ~isreal(D)|i==0
  warning('ODESPX: x is not saddle point or stable steady-state.')
  x = [];
  return
end

[~,i] = min(diag(D));
delx = 0.0001*V(:,i);
t = nodeunif(n,0,-T);
h = [0;diff(t)];

xsp = zeros(n,2);
x = x0+delx;
xsp(1,:) = x;
for i = 2:n
  hh = h(i);
  f1 = feval(f,x)*(hh/2);
  f2 = feval(f,x+f1)*hh;
  f3 = feval(f,x+f2/2)*hh;
  f4 = feval(f,x+f3)*(hh/2);
  x = x+(f1+f2+f3+f4)/3;
  xsp(i,:,:) = x;
end
xsp(i+1:n,:) = [];
xsp = real(xsp);

xsn = zeros(n,2);
x = x0-delx;
xsn(1,:) = x;
for i = 2:n
  hh = h(i);
  f1 = feval(f,x)*(hh/2);
  f2 = feval(f,x+f1)*hh;
  f3 = feval(f,x+f2/2)*hh;
  f4 = feval(f,x+f3)*(hh/2);
  x = x+(f1+f2+f3+f4)/3;
  xsn(i,:,:) = x;
end
xsn(i+1:n,:) = [];
xsn = real(xsn);

x = [xsn(end:-1:1,:);x0';xsp];
i = all(isfinite(x),2);
x = x(i,:);