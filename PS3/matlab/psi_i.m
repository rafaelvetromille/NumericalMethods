function psi = psi_i(i, k, nk)

alpha = 1/3;
beta  = 0.987;
delta = 0.012;

k_ss  = ((1/beta + delta - 1)*(1/alpha))^(1/(alpha - 1));
kgrid = linspace(0.75*k_ss, 1.25*k_ss, nk);

if i == 1
    
    if k >= kgrid(i) && k <= kgrid(i+1)   
        psi = (kgrid(i+1) - k)/(kgrid(i+1) - kgrid(i));   
    else 
        psi = 0;
    end
    
elseif i == nk
    
    if k <= kgrid(i) && k >= kgrid(i-1)     
        psi = (k - kgrid(i-1))/(kgrid(i) - kgrid(i-1));
    else    
        psi = 0;
    end
    
else
    
    if k >= kgrid(i-1) && k <= kgrid(i)
        psi = (k - kgrid(i-1))/(kgrid(i) - kgrid(i-1));  
    elseif k >= kgrid(i) && k <= kgrid(i+1)
        psi = (kgrid(i+1) - k)/(kgrid(i+1) - kgrid(i)); 
    else
        psi = 0;
    end
end

end

