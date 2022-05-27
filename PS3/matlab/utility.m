function u = utility(c,mu)

if mu==1
    
    u = log(c);
    
else
    
    u = (1/(1-mu))*(c.^(1-mu) - 1);

end