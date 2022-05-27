function [Tv, idx] = iteration_multigrid(v, multigrid, kgrid_mg1, func, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

nmg       = length(multigrid);
nz        = length(zgrid);
kss       = (1/alpha*(1/beta + delta - 1))^(1/(alpha-1));            % Valor do capital de estado estacionário;

for i = 1:nmg
           
    [Tv, idx] = func(v, kgrid_mg1, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
    
    if i < nmg
        % Definimos o primeiro grid (mais fino) do processo de multigrid
        nki       = multigrid(i+1); % número de elementos dos próximos grids da iteração
        kgrid_mgi = linspace(0.75*kss, 1.25*kss, nki);
        v         = zeros(nki,nz);
        
        for iz = 1:nz
            %Agora, interpolamos o kmg1 (de menor dimensão) através de Tv num
            %vetor de dimensão nki > nk1            
            v(:,iz) = interp1(kgrid_mg1, Tv(:,iz), kgrid_mgi, 'linear');   % linear or spline?
        end
        
        % Substituimos pelo outro grid e repetimos o processo
        kgrid_mg1 = kgrid_mgi;
    end    
end
end