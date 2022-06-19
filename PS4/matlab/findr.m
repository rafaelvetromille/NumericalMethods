function z = findr(r,parameters,max_iter,penalty,tol)

beta  = parameters.beta;
mu    = parameters.mu;

agrid = parameters.agrid; na = numel(agrid);
zgrid = parameters.zgrid; nz = numel(zgrid);

a     = parameters.a; ns = numel(a);
z     = parameters.z;
P     = parameters.P;

%%%%% reward function

c    = zeros(ns,na);

for j=1:na
    
    c(:,j) = (1+r)*a + exp(z) - agrid(j);
  
end

%%%%% check consumption non-negative

violations = (c<0);

c = c.*(c>=0) + eps;  % soma um valor muito pequeno para não ficar zero

% UTILITY FUNCTION

if mu==1
    u  = @(c)log(c);
    u0 = log(c) - penalty*violations;
else    
    u  = @(c)(c.^(1-mu)-1)./(1-mu);
    u0 = (1/(1-mu))*(c.^(1-mu) - 1) - penalty*violations;
end

% INITIALIZE VALUE FUNCTION

v1        = log(0.5*z)/(1-beta);                       % initial guess 1
%v2       = (1/(1-mu))*((r.*a+y).^(1-mu) - 1);         % initial guess 2
%v3       = reshape(Vguess,[],1);                      % initial guess 3
v4       = nodeunif(1809, 0, 0);                       % initial guess 4

% ITERATION OF VALUE FUNCTION

v       = v4;  % choose a initial guess (all of them give us the same result)

for i=1:max_iter

RHS      = u0+beta*kron(P,ones(na,1))*reshape(v,na,nz)';

[Tv,argmax] = max(RHS,[],2);

%%%%% policy that attains the maximum

g = a(argmax);

%%%%% check if converged

error = max(abs(Tv-v));

%fprintf('%4i %6.2e \n',[i, error]);

if norm(Tv-v,inf)<tol, break, end

%%%%% if not converged, update guess

v = Tv;

end

%%%%% market clearing condition

%%%%% construct transition matrix for the state s=(a,y)

A = zeros(ns,na);
Q = zeros(ns,ns);

PP = kron(P,ones(na,1));

for s=1:ns
    
    A(s,:) = (agrid==g(s))';  %% puts a 1 if g(s)=a 
    
    Q(s,:) = kron(PP(s,:),A(s,:));
    
end

%%%%% compute stationary distribution

[eig_vectors,eig_values] = eig(Q'); 
[~,arg] = min(abs(diag(eig_values)-1)); 
unit_eig_vector = eig_vectors(:,arg); 

lambda = unit_eig_vector/sum(unit_eig_vector); 

%%%%% check market clearing

z = sum(lambda.*g);

fprintf('[r, z=F(r)] %2.9f  %3.8f  \n', [r, z]);

