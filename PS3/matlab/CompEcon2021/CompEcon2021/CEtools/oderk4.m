%% ODERK4
%
%  Solves initial value problem using fourth-order Runge-Kutta method
%
%  Specifically, solves first-order ordinary differential equation
%    x'(t) = f(x(t)), t in [0,T]  s.t. x(0) = x0
%  Here, x is d.1 vector-valued function defined on time domain [0,T] and
%  x' is its d.1 vector-valued derivative with respect to time.
%
%  Usage
%    [x,t] = oderk4(f,x0,T,n,varargin)
%  Let
%    d  = dimension of state process x
%    k  = number of initial values passed to odekr4
%  Input
%    f         : velocity function
%    x0        : d.k initial values
%    T         ; time horizon
%    n         : number of time nodes
%    varargin  : optional parameters passed to f
%  Output
%    x         : n.d.k values at time nodes
%    t         : n.1 equally spaced time nodes on [0,T]

%  Velocity Function
%    User-supplied function that returns velocity at states according to
%    Format
%      v = f(x,varargin)
%    Input
%      x         : d.k states
%      varargin  : optional parameters passed to f
%    Output
%      v         : d.k velocities

%  Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,t] = oderk4(f,x0,T,n,varargin)

if nargin<4 || isempty(n)
  n = 1000;
end

t = nodeunif(n,0,T);
[d,k] = size(x0);
n = length(t);
x = zeros(n,d,k);
x(1,:,:) = x0;
h = [0;diff(t)];
for i = 2:n
  hh = h(i);
  f1 = feval(f,x0,varargin{:})*(hh/2);
  f2 = feval(f,x0+f1,varargin{:})*hh;
  f3 = feval(f,x0+f2/2,varargin{:})*hh;
  f4 = feval(f,x0+f3,varargin{:})*(hh/2);
  x0 = x0+(f1+f2+f3+f4)/3;
  x(i,:,:) = x0;
end
x = squeeze(real(x));