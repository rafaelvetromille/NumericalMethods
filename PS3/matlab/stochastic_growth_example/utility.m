function u = utility(c,sigma)

if sigma==1,
    
    u = log(c);
    
else
    
    u = (1/(1-sigma))*(c.^(1-sigma) - 1);

end