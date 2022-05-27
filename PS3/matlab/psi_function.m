function psi = psi_function(i, k, parameters)

% No i = 1, usando a segunda linha, pois não tenho i-1.
if i == 1
    
    if k >= parameters.kgrid(i) && k <= parameters.kgrid(i+1)
        psi = (parameters.kgrid(i+1) - k)/(parameters.kgrid(i+1) - parameters.kgrid(i));
    else
        psi = 0;
    end
    
% No i = nz, usando a primeira linha, pois não tenho i+1.
elseif i == length(parameters.kgrid)
    
    if k >= parameters.kgrid(i-1) && k <= parameters.kgrid(i)
        psi = (k - parameters.kgrid(i-1))/(parameters.kgrid(i) - parameters.kgrid(i-1));
    else
        psi = 0;
    end
    
% Nos demais segue-se a regra natural descrita na função.
else 
    
    if k >= parameters.kgrid(i-1) && k <= parameters.kgrid(i)
        psi = (k - parameters.kgrid(i-1))/(parameters.kgrid(i) - parameters.kgrid(i-1));
    elseif k >= parameters.kgrid(i) && k <= parameters.kgrid(i+1)
        psi = (parameters.kgrid(i+1) - k)/(parameters.kgrid(i+1) - parameters.kgrid(i));
    else
        psi = 0;
    end
    
end
end

