function psi = psi_function(i, k, n_elements, parameters)

kmin      = parameters.kgrid(1); 
kmax      = parameters.kgrid(end);

kgrid_new = linspace(kmin, kmax, n_elements);

% No i = 1, usando a segunda linha, pois não tenho i-1.
if i == 1
    
    if k >= kgrid_new(i) && k <= kgrid_new(i+1)
        psi = (kgrid_new(i+1) - k)/(kgrid_new(i+1) - kgrid_new(i));
    else
        psi = 0;
    end
    
% No i = nz, usando a primeira linha, pois não tenho i+1.
elseif i == length(kgrid_new)
    
    if k >= kgrid_new(i-1) && k <= kgrid_new(i)
        psi = (k - kgrid_new(i-1))/(kgrid_new(i) - kgrid_new(i-1));
    else
        psi = 0;
    end
    
% Nos demais segue-se a regra natural descrita na função.
else 
    
    if k >= kgrid_new(i-1) && k <= kgrid_new(i)
        psi = (k - kgrid_new(i-1))/(kgrid_new(i) - kgrid_new(i-1));
    elseif k >= kgrid_new(i) && k <= kgrid_new(i+1)
        psi = (kgrid_new(i+1) - k)/(kgrid_new(i+1) - kgrid_new(i));
    else
        psi = 0;
    end
    
end
end

