%% LQAPPROX2

% [sstar,xstar,pstar,vstar,X,P,G,xlq,plq,vlq] = lqapprox2(delta,f,g,s0,x0,s,params)
function [sstar,xstar] = lqapprox2(delta,f,g,s0,x0)

maxit = 100;

% ns = size(s,1);
ds  = length(s0);

sstar = s0;
xstar = x0;

for it=1:maxit
  % Compute LQ Approximation Around Steady-State Using lqsolve
  sstarold = sstar;
  xstarold = xstar;
  f0  = f([s0;x0]);
  fJ  = fdjac(f,[s0;x0]);
  fH  = fdhess(f,[s0;x0]);
  fs  = fJ(1:ds);
  fx  = fJ(ds+1:end);
  fss = fH(1:ds,1:ds);
  fsx = fH(1:ds,ds+1:end);
  fxs = fH(ds+1:end,1:ds);
  fxx = fH(ds+1:end,ds+1:end);
  gJ  = fdjac(g,[s0;x0]);
  gs  = gJ(:,1:ds);
  gx  = gJ(:,ds+1:end);
  f0  = f0 - fs*s0 - fx*x0 + 0.5*s0'*fss*s0 + s0'*fsx*x0 + 0.5*x0'*fxx*x0;
  fs  = fs - s0'*fss - x0'*fxs;
  fx  = fx - s0'*fsx - x0'*fxx;
  s0  = s0 - gs*s0 - gx*x0;
  [~,~,~,sstar,xstar] = lqsolve(f0,fs,fx,fss,fsx,fxx,s0,gs,gx,delta);
  change = norm(sstar-sstarold)+norm(xstar-xstarold);
  fprintf('\n %2i %9.3e \n',it,change)
  if change<1.e-8
    break
  else
    s0 = sstar;
    x0 = xstar;
  end
end

% % Get derivatives
% [f0,fx,fxx] = feval(func,'f',s0,x0,1,1,estar,params{:});
% [g0,gx]     = feval(func,'g',s0,x0,1,1,estar,params{:});
% fs  = fjac(func,2,'f',s0,x0,1,1,estar,params{:});
% fxs = fjac(func,[2,2],'f',s0,x0,1,1,estar,params{:});
% fss = fhess(func,2,'f',s0,x0,1,1,estar,params{:});
% gs  = fjac(func,2,'g',s0,x0,1,1,estar,params{:});
% 
% % Reshape to ensure conformability
% s0 = s0';
% x0 = x0';
% fs  = reshape(fs , 1,ds);
% fx  = reshape(fx , 1,dx);
% fss = reshape(fss,ds,ds);
% fxs = reshape(fxs,dx,ds);
% fxx = reshape(fxx,dx,dx);
% g0  = reshape(g0 ,ds, 1);
% gx  = reshape(gx ,ds,dx);
% gs  = reshape(gs ,ds,ds);
% fsx = fxs';
% f0 = f0 - fs*s0 - fx*x0 + 0.5*s0'*fss*s0 + s0'*fsx*x0 + 0.5*x0'*fxx*x0;
% fs = fs - s0'*fss - x0'*fxs;
% fx = fx - s0'*fsx - x0'*fxx;
% g0 = g0 - gs*s0 - gx*x0;
% 
% % Solve Ricatti equation using QZ decomposition
% A = [eye(ds)          zeros(ds,dx+ds);
%      zeros(dx,ds+dx) -delta*gx'      ;
%      zeros(ds,ds+dx)  delta*gs'     ];
% B = [ gs   gx  zeros(ds,ds);
%      fsx' fxx  zeros(dx,ds);
%     -fss -fsx  eye(ds)     ];
% [~,~,~,Z] = qzordered(A,B);
% C  = real(Z(ds+1:end,1:ds)/Z(1:ds,1:ds));
% X  = C(1:dx,:);
% P  = C(dx+1:end,:);
% 
% % Compute steady-state state, action, and shadow price
% t = [fsx' fxx delta*gx';fss fsx delta*gs'-eye(ds);gs-eye(ds) gx zeros(ds,ds)]\[-fx';-fs'; -g0];
% sstar = t(1:ds);
% xstar = t(ds+1:ds+dx);
% pstar = t(ds+dx+1:ds+dx+ds);
% vstar = (f0+fs*sstar+ fx*xstar+0.5*sstar'*fss*sstar+sstar'*fsx*xstar+0.5*xstar'*fxx*xstar)/(1-delta);

% % Compute lq-approximation optimal policy and shadow price functions at
% % state nodes
% sstar = sstar';
% xstar = xstar';
% pstar = pstar';
% s = s-sstar(ones(ns,1),:);
% xlq = xstar(ones(1,ns),:) + s*X';
% plq = pstar(ones(1,ns),:) + s*P';
% vlq = vstar + s*pstar' + 0.5*sum(s.*(s*P'),2); 



%
%  Derives linear-quadratic approximation to the solution of the
%  determministic infinite-horizon Bellman equation 
%
%    V(s) = max_x f(s,x) + delta*V(g(s,x)).

%  where s is the state variable, x is the action variable, f is the reward
%  function, and t is transition function g.
%
%  Specifically, computes the steady-state sstar,xstar, pstar, and vstar;
%  replaces f and g, respectively, with their Taylor quadratic and linear
%  approximations around the steady state;, and solves the resulting
%  linear-quadratic control problem to generate the approximations
%
%    Optimal policy function
%       x(s) = xstar + X*(s-sstar)
%    Shadow price function
%       p(s) = pstar + P*(s-sstar)
%    Value function
%       V(s) = vstar + pstar*(s-sstar) + 0.5*(s-sstar)'*P*(s-sstar)
%    Controlled state transition function
%      snext = sstar + G*(s-sstar).
%
%  Usage
%    [sstar,xstar,pstar,vstar,X,P,G,xlq,plq,vlq] = lqapprox2(delta,f,g,s0,x0,s,params)
%  Let
%    ns = number of state nodes input
%    ds = dimension of state s
%    dx = dimension of action x
%  Input
%    delta  : discount factor
%    f      : reward function
%    g      : transition function
%    s0     : 1.ds guess for the steady-state state
%    x0     : 1.dx guess for the steady-state action
%    s      : ns.ds state nodes for approximations
%    params : user-supplied list of parameters passed to f and g
%  Output
%    sstar  : 1.ds  steady-state state
%    xstar  : 1.dx  steady-state action
%    pstar  : 1.ds  steady-state shadow price
%    vstar  : 1.ds  steady-state value
%    X      : dx.ds slope of policy function
%    P      : ds.ds slope of shadow price function
%    G      : ds.ds sloppe of controlled state transition function
%    xlq    : ns.dx optimal lq approximation actions at state nodes s
%    plq    : ns.ds optimal lq approximation shadow prices at state nodes s
%    vlq    : ns.1  optimal lq approximation values at state nodes s
% 
%  Function Files
%    User-supplied functions f and g take the form
%      fval = f(s,x,params)
%      gval = g(s,x,params)
%  where
%      s         : 1.ds vector of states
%      x         : 1.dx vector of actions
%      params    : user-supplied list of function parameters
%      fval      : 1.1 reward function value
%      gval      : 1.ds transition function values

%  Copyright(c) 2021
%   Mario J. Miranda - miranda.4@osu.edu