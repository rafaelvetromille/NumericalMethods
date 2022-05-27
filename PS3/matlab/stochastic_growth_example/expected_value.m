function y = expected_value(c,s,parameters,a,fspace)

alpha = parameters.alpha;
delta = parameters.delta;

rhoz  = parameters.rhoz;
ez    = parameters.ez;
wz    = parameters.wz;

k = s(:,1); 
z = s(:,2); zmin = z(1); zmax = z(end);

Ev = 0;

kprime = z.*(k.^alpha) + (1-delta)*k-c;

for j=1:numel(ez),

zprime = max(min(z.^rhoz.*exp(ez(j)), zmax), zmin);

sprime = [kprime,zprime];

Ev = Ev+wz(j)*funeval(a,fspace,sprime);

end

y = Ev;
