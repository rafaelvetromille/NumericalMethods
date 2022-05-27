function [Tv, idx] = iteration_monotonicity_and_concacity(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

iter  = 0;
error = 1;

nk    = length(kgrid);
nz    = length(zgrid);

Tv    = zeros(nk,nz);
idx   = zeros(nk,nz);

% Itera��o da Fun��o valor - Concavidade c/ Acelerador
% ----------------------------------------------------
tic;
while  error > tol && iter <= max_iter
    
    for iz = 1:nz
        
        z = zgrid(iz);
        
        % o 'start' vai sendo atualizado a cada itera��o
        % come�amos pelo caso igual ao anterior
        % ---------------------------------------------------
        start = 1;
        
        for ik = 1:nk
            
            k = kgrid(ik);
            c = zeros(nk,1);
            H = zeros(nk,1);
            
            %%% Primeiro, computar para ikprime = 1
            ikprime = start;
            
            %%% Computar o consumo em ikprime = 2
            c(ikprime) = z*(k^alpha) + (1-delta)*k - kgrid(ikprime);
            
            %%% Verificar consumos negativos
            if c(ikprime) >= 0
                u0 = (c(ikprime)^(1-mu)-1)/(1-mu);
            else
                u0 = -Inf;
            end
            
            %%% Calcular a esperan�a
            Ev = 0;
            for jz = 1:nz
                Ev = Ev + P(iz,jz)*v(ikprime,jz);
            end
            
            %%% Computar a fun��o valor
            H(ikprime) = u0 + beta.*Ev;
            
            %%% Agora, computar ikprime de 2 em diante
            for ikprime = (start+1):nk
                
                %%% Computar o consumo a partir de ikprime = 2
                c(ikprime) = zgrid(iz)*(k^alpha) + (1-delta)*k - kgrid(ikprime);
                
                %%% Verificar consumos negativos
                if c(ikprime) >= 0
                    u0 = (c(ikprime)^(1-mu)-1)/(1-mu);
                else
                    u0 = -Inf;
                end
                
                %%% Calcular a esperan�a (bem mais r�pido)
                Ev = 0;
                for jz = 1:nz
                    Ev = Ev + P(iz,jz)*v(ikprime,jz);
                end
                
                %%% Computar a fun��o valor
                H(ikprime) = u0 + beta.*Ev;
                
                % Calcular a convidade aqui! Veja que n�o h� uso do max.
                if H(ikprime) < H(ikprime-1)
                    break
                end
            end
            
            %%% Mudan�a caso nunca tenha acionado o break (ou seja, o m�ximo � no �ltimo t mesmo)
            if  H(ikprime) <= H(ikprime-1)
                
                %%% Atualizar a fun��o valor e a matriz de �ndices
                Tv(ik,iz)  = H(ikprime-1); idx(ik,iz) = ikprime-1;
                
            elseif H(ikprime) >= H(ikprime-1)
                
                %%% Atualizar a fun��o valor e a matriz de �ndices
                Tv(ik,iz)  = H(ikprime); idx(ik,iz) = ikprime;
                
            end
            
            %%% O novo 'start' ser� a partir da decis�o da itera��o anterior
            start = idx(ik,iz);
            
        end
    end
    
    %%% Computar o erro
    error = max(max(abs(Tv - v)));
    iter = iter + 1;                    % atualizamos o contador de passos da itera��o
    v = Tv;                             % atualizamos o chute de V pelo TV encontrado
    
    %%% Printar o n�mero de itera��es e o erro
    fprintf('Error %4i %6.2e \n', [iter, error]);
    pause(0.001);
    clc;

end
time = toc;

% Imprimir os resultados (tempo decorrido, n� de itera��es e error)
% -----------------------------------------------------------------
last_iter = iter;

if error<tol && last_iter<=max_iter
    
    fprintf('\n');
    fprintf('Method 4. Iteration of Value Function (Monotonicity & Concavity) \n\n');
    fprintf('Time (in seconds)      = %3.03f \n',time);
    fprintf('Number of iterations   = %3.00f \n',last_iter);
    fprintf('Error                  = %3.10f \n',error);
    fprintf('\n');
    
end
end