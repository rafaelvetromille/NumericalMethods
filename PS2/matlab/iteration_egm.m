function C1 = iteration_egm(kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

iter   = 0; 
error  = 1; 

nk     = length(kgrid);
nz     = length(zgrid);

C0     = zeros(nk,nz);
C1     = zeros(nk,nz);

%k_end  = zeros(nk,nz);  % grid endógeno do k
k_exo  = kgrid;          % grid exógeno do k

k_end           = zeros(nk,1);
k_prime_end     = zeros(nk,1);
           
% chute inicial pra o consumo
for iz=1:nz
    for ik=1:nk
        C0(ik,iz) = zgrid(iz)*(k_exo(ik)^alpha) + (1-delta)*k_exo(ik) - k_exo(ik);
    end
end

tic;
while  error > tol && iter <= max_iter
    for iz = 1:nz
        z = zgrid(iz);
        for ik = 1:nk    
            % calcular a esperança (bem mais rápido) - Chris Edmond
            Ec = 0;
            for w = 1:nz
                total = (C0(ik,w)^(-mu))*(zgrid(w)*alpha*(k_exo(ik)^(alpha-1)) + (1-delta));
                Ec    = Ec + beta*P(iz,w)*total;
            end
            
            % Lado direito da equação de Euler
            RHS = k_exo(ik) + Ec^(-1/mu);
            
            % pela CPO, deveríamos ter que c = aux3
            % pra cada z,k' vamos procurar k que zera a função
            % pelo método da bisecção
            f               = @(k)(z*(k^alpha) + (1-delta)*k - RHS);            
            k_end(ik)       = bisection_method(f, 0, 10e5);
            k_prime_end(ik) = k_exo(ik);
            
        end
        
        % interpolo para achar k' ótimo no grid exógeno correspondente ao grid
        % endógeno
        k_prime_exo = interp1(k_end, k_exo, k_exo, 'spline', 'extrap');  % linear
        %k_prime_exo = spline(k_end, k_prime_end, k_exo);
        
        % calcular a função política do consumo (não pode ser negativo)
        for jk = 1:nk
            C1(jk,iz) = z*(k_exo(jk)^alpha) + (1-delta)*k_exo(jk) - k_prime_exo(jk);
        end
        
    end
    
    % Computar o erro
    error = max(max(abs(C1-C0)));
    iter  = iter+1;
    
    % atualizamos o chute
    C0    = C1;
    
    % printar o número de iterações e o erro
    fprintf('Error %4i %6.2e \n', [iter, error]);
    
end
time = toc;

% Imprimir os resultados (tempo decorrido, nº de iterações e error)
% -----------------------------------------------------------------
last_iter = iter;

if error<tol && last_iter<=max_iter
    
    fprintf('\n');
    fprintf('Method 10. Iteration of Value Function (Endogenous Grid) \n\n');
    fprintf('Time (in seconds)      = %3.03f \n',time);
    fprintf('Number of iterations   = %3.00f \n',last_iter);
    fprintf('Error                  = %3.10f \n',error);
    fprintf('\n');
    
end