function [Tv, idx] = iteration_bruteforce(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol)

iter  = 0;
error = 1;

nk    = length(kgrid);
nz    = length(zgrid);

Tv    = zeros(nk,nz);
idx   = zeros(nk,nz);

tic;
while (error > tol && iter <= max_iter)
    for iz = 1:nz
        
        z = zgrid(iz);
        
        for ik = 1:nk
            
            % Pegando o capital ik e criando um vetor de consumo
            % -----------------------------------------------------------
            k = kgrid(ik);                  % por algum motivo aumenta a performance do c�gido
            c = zeros(nk,1);
            H = zeros(nk,1);
            
            % Esse loop define o vetor de consumo de forma a garantir que
            % nenhuma entrada ser� negativa: a cada elemento, tiramos o
            % m�ximo entre o consumo e zero
            % -----------------------------------------------------------
            for ikprime = 1:nk
                
                %%% Computar o consumo
                c(ikprime)  = z.*(k^alpha) + (1-delta).*k - kgrid(ikprime);
                
                %%% Verificar consumos negativos
                if c(ikprime) >= 0
                    u0 = (c(ikprime)^(1-mu)-1)/(1-mu);
                else
                    u0 = -Inf;
                end
                
                %%% Computar a esperan�a (Chris Edmond)
                Ev = 0;
                for jz = 1:nz
                    Ev = Ev + P(iz,jz)*v(ikprime,jz);
                end
                
                %%% Computar a fun��o valor
                H(ikprime) = u0 + beta.*Ev;
                
            end
            
            %%% Computar a fun�ao valor (matriz)
            [Tv(ik,iz), idx(ik,iz)] = max(H);
            
        end
    end
    
    % Uma vez finalizado esse loop, iremos checar se Tv est� suficientemente
    % perto de v. Como � imposs�vel que Tv = v, usamos uma toler�ncia
    % como uma aproxima��o razo�vel.
    
    % Note que:
    % abs(TV-V): � uma matriz 500x7 composta pelo valor absoluto da diferen�a das duas matrizes elemento a elemento
    % max(abs(TV-V)): retorna o maior elemento de cada coluna, entao � 1x7
    % max(max(abs(TV-V))): retorna o maior elemento desse vetor (1x1)
    error = max(max(abs(Tv - v)));
    iter  = iter + 1;                   % atualiza a  itera��o
    v     = Tv;                         % atualiza o chute de V pelo novo TV encontrado
    
    % Imprime a itera��o e o erro
    fprintf('Error %4i %6.2e \n',[iter, error]);
    pause(0.01)
    clc
end
time = toc;

% Imprimir os resultados (tempo decorrido, n� de itera��es e error)
% -----------------------------------------------------------------
last_iter = iter;

if error<tol && last_iter<=max_iter
    
    fprintf('\n');
    fprintf('Method 1. Iteration of Value Function (Brute Force) \n\n');
    fprintf('Time (in seconds)      = %3.03f \n',time);
    fprintf('Number of iterations   = %3.00f \n',last_iter);
    fprintf('Error                  = %3.10f \n',error);
    fprintf('\n');
    
end
end