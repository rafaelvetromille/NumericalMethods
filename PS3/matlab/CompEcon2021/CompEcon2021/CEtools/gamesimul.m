%% GAMESIMUL
%
%  Simulates discrete time two person dynamic game
%
%  Simulates the Markov process followed by the optimized state variable of
%  a discrete-time continuous-state continuous-action infinite-horizon
%  multi-agent dynamic game model solved using gamesolve.
%
%  Usage
%    [ssim,xsim] = gamesimul(model,ss,nper,sres,xres,esim)
%  Input
%    model   : model structure variable
%    ss      : k by d vector of initial states
%    nper    : number of simulated time periods
%    sres    : coordinates of the evaluation grid (from dpsolve)
%    xres    : optimal control function values at grid defined by sres
%    esim    : nrep.nper.de simulated shocks
%  Output
%    ssim    : k by d by nper+1 vector of simulated states       
%    xsim    : k by d by nper+1 vector of simulated actions
%    For finite horizon problems, xsim is k by d by nper
%    If d=1, ssim and xsim are k by nper+1
%  See
%    gamesolve

%  Copyright(c) 1997-2021
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function [ssim,xsim] = gamesimul(model,ss,nper,sres,xres,esim)
func   = model.func;
params = model.params;
nrep   = size(ss,1);
ds     = size(ss,2);
st     = gridmake(sres);
dx     = ds*length(xres(:))/length(st(:));
ssim   = zeros(nrep,nper+1,ds);
xsim   = zeros(nrep,nper+1,dx);

nx = numel(xres)/dx;
xres = reshape(xres,nx,dx);
for ip=1:nper+1
  xx = minterp(sres,xres,ss);
  ssim(:,ip,:) = ss;
  xsim(:,ip,:) = xx;
  ss = feval(func,'g',1,ss,xx,squeeze(esim(:,ip,:)),params{:});
end

if ds==1; ssim=squeeze(ssim); end
if dx==1; xsim=squeeze(xsim); end