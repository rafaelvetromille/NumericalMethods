function x = solve_brent(f, k,parameters,coefficients,fspace,a, b, tol)

v = 1/2*(a + b);       % initial guess

alpha1 = (3-sqrt(5))/2;

w    = v;                                                                     % key three point in parabolic steps 
x    = v;   
fx   = -feval(f, x, k,parameters,coefficients,fspace);

fv   = fx;
fw   = fx;
d    = zeros(size(a));                                                        % just put a number 
e    = zeros(size(a));                                                        % choosen s.t. jumps directly to GS 
term = 1;

while any(term)

   xm   = 0.5*(a+b);                                                          % middle point            
   tol1 = tol*abs(x) + eps;
   tol2 = 2*tol1;
   term = abs(x-xm) > tol2 - 0.5*(b-a);                                       % if any term is bigger than 1 the optimizer continue
   
   % parabolic step
   
   par = abs(e) > tol1;                                                       % First control to parabolic, step has to be big enough     
   m   = (x-w).*(fx-fv);
   q   = (x-v).*(fx-fw);
   p   = (x-v).*q-(x-w).*m;                                                   % numerator in the parabolit fit
   q   = 2.*(q-m);                                                            % denominator in the parabolit fit      
   
   p   = p.*(-(q > zeros(size(q)))*2+1);                                      % if q is positive make sure p is negative for the x-d step   
   q   = abs(q);                                                              % IMPORTANT: SIGNS FOR THE PARABOLIC STEPS 
   
   etemp = e;                                                              
   e_p   = d;                                                                
   
   par  = par & ~((abs(p)>=abs(0.5.*q.*etemp))|(p<=q.*(a-x))|(p>= q.*(b-x))); %  If this step fails then you jump to golden
   
   d = p./(q+(q==0));                                                        % never divide by zero, if this is the case, par is 0 by the first term 
   u = x + d;
   aux = u - a < tol2 | b - u < tol2;
   d_p = aux.*sign(xm-x).*tol1+(~aux).*d;
   
   % golden search step 
   
   e_g = (a-x).*(x >= xm)+(b-x).*(x<xm);                                     % golden search in choosing d
   d_g = alpha1*e_g;                                                         % goldent step
   
   % choose the step (d) and e of the parabolic (if it past all the check, if not golden)
   d = par.*d_p+(~par).*d_g;
   e = par.*e_p+(~par).*e_g;
   
   % this steps are independent of parabolic or golden steps 
   
   u  = (x+d).*(abs(d)>=tol1)+(abs(d)<tol1).*(x+sign(d).*tol1);
   fu = -feval(f, u, k,parameters,coefficients,fspace);
   
   % redefine the upper and lower limit in the optimization: Both for
   % golden and for parabolic 
   
   aux1 = fu <= fx;   
   aux2 = u  >=  x;
      
   a = aux1.*(aux2.*x+(~aux2).*a)+(~aux1).*((~aux2).*u+aux2.*a);
   b = aux1.*((~aux2).*x+aux2.*b)+(~aux1).*(aux2.*u+(~aux2).*b);
   
   % adjusting the three point in the parabolic step as a function fu & fx
   
   aux3 = (fu<=fw)|(w==x);
   aux4 = (fu<=fv)|(x==v)|(v==w);
   
   v  = aux1.* w+(~aux1).*(aux3.*w +(~aux3).*((aux4).* u+(~aux4).* v));     
   fv = aux1.*fw+(~aux1).*(aux3.*fw+(~aux3).*((aux4).*fu+(~aux4).*fv));
   
   w  = aux1.* x+(~aux1).*(aux3.* u+(~aux3).*w);      
   fw = aux1.*fx+(~aux1).*(aux3.*fu+(~aux3).*fw);
   
   x  = aux1.*u +(~aux1).*x;
   fx = aux1.*fu+(~aux1).*fx;
   
end

