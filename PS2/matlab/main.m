%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  M�todos Num�ricos - Lista 2                            %
%                            Programa escrito por: Rafael Vetromille                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In�cio formal do documento
% --------------------------
clearvars;
close all;
clc;

%%

% Exercise 3
% ----------

% De agora em diante, use o modelo completo com incerteza. Resolva o
% problema no computador utilizando o m�todo da itera��o da fun��o
% valor padr�o. Para tanto, voc� ter� que discretizar suas vari�veis de
% estado. Para a choque de TFP, utilize o m�todo de Tauchen (1986) com 7 pontos.

% Calibra��o
% ----------
beta     = 0.987;                        % Fator de desconto; 
mu       = 2;                            % Coeficiente de avers�o relativa ao risco;
alpha    = 1/3;                          % Sahre do capital na fun��o de produ��o; 
delta    = 0.012;                        % Taxa de deprecia��o;        
const    = 0;                            % Intercept term of the AR(1) process
rho      = 0.95;                         % Par�metro de persist�ncia dos choques de renda;
sigma    = 0.007;                        % Vari�ncia do erro do processo de renda.

nz       = 7;                            % N�mero de estados ex�genos para a vers�o discreta do processo de renda AR(1); 
m        = 3;                            % N�mero de desvios em rela��o � m�dia do processo AR para definir a amplitude do espa�o de estados. 

% Capital de estado estacion�rio
% ------------------------------
kss      = (1/alpha*(1/beta + delta - 1))^(1/(alpha-1));            % Valor do capital de estado estacion�rio;
fprintf('Com a calibra��o proposta, o capital de estado estacion�rio (kss) � %3.2f. \n\n', kss);

% M�todo de Tauchen
% -----------------

% FUNCTION: [s, Pi] = mytauchen(const, m, rho, sig, N)
%
% INPUTS:
%   const:  escalar, intercept term of the AR(1) process
%   rho:    escalar, coeficiente do AR(1)
%   sig:    escalar, desvio padr�o dos choques de TFP
%   N:      escalar, n�mero de pontos para o processo discretizado
%
% OUTPUTS:
%   S:      column vector of size Nx1, contains all possible states in ascending order
%   Pi:     matrix of size NxN, contains the transition proabilities. Rows are current state and columns future state
% -------------------------------------------------------------------------------------------------------------------
[S, P]   = mytauchen(const, m, rho, sigma, nz);

% Para o grid de capital, use 500 pontos linearmente espa�ados
% no intervalo [0.75*kss, 1.25*kss]. Eu recomendo fortemente que voc� n�o
% use o m�todo da "for�a-bruta" para encontrar a fun��o pol�tica.
% -----------------------------------------------------------------------
zmin     = -m*sqrt(sigma^2/(1 - rho^2));
zmax     = +m*sqrt(sigma^2/(1 - rho^2));

nk       = 500; 
kmax     = 1.25*kss;
kmin     = 0.75*kss;

kgrid    = linspace(kmin, kmax, nk)';
zgrid    = exp(linspace(zmin, zmax, nz))';

% Par�metros num�ricos
% --------------------
tol      = 10e-6;                   % Toler�ncia do erro.
iter     = 0;                       % N�mero de itera��es;
error    = 1;                       % Primeiro error que dar� o start no while;
max_iter = 10000;                   % N�mero m�ximo de itera��es;

% Fun��o de utilidade (caso geral)
% --------------------------------
if mu == 1
    u  = @(c)log(c);
else    
    u  = @(c)(c.^(1-mu)-1)./(1-mu);
end

% Guess Inicial (Benjamin Moll) - poderia ser zero, mais � menos eficiente
v = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        v(ik,iz) = u(zgrid(iz).*kgrid(ik).^alpha + (1-delta).*kgrid(ik) - kgrid(ik))./(1-beta);
    end
end

%% Brute Force
clc;

% M�todo 1: Itera��o da Fun��o Valor pela For�a Bruta
tic;
[v1, idx1] = iteration_bruteforce(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(1) = toc;

%% Fun��es pol�ticas

% Extrair as fun��es pol�tica k' = g(k,z):
g1 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        g1(ik,iz) = kgrid(idx1(ik,iz));
    end
end

% Extrair as fun��es pol�tica c = c(k,z):
c1 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        c1(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha)+(1-delta)*kgrid(ik)-g1(ik,iz);
    end
end

%% Gr�fico - Fun��o Valor (2D & 3D) - Brute Force
figure(1)
plot_value_function(v1, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Capital) (2D & 3D) - Brute Force
figure(2)
plot_capital_policy_function(g1, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Consumo) (2D & 3D) - Brute Force
figure(3)
plot_consumption_policy_function(c1, kgrid, zgrid);

%% Euler Equation Erros - Brute Force
figure(4)
EEE1 = euler_equation_erros(idx1, c1, kgrid, zgrid, P, alpha, beta, delta, mu);

% Verificar o mpinimo valor do EEE
max(max(EEE1))    % maior valor
min(min(EEE1))    % menor valor

%% Brute Force w/ Accelerator
clc;

% M�todo 1: Itera��o da Fun��o Valor pela For�a Bruta c/ Acelerador 
tic;
[v2, idx2] = iteration_bruteforce_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(2) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx2);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Brute Force com Acelerador s�o %s. \n\n', teste);

%% Monotonicity 
clc;

% M�todo 2: Itera��o da Fun��o Valor pela Monotonicidade
tic;
[v3, idx3] = iteration_monotonicity(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(3) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx3);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade s�o %s. \n\n', teste);

%% Monotonicity w/ Accelerator
clc;

% M�todo 2: Itera��o da Fun��o Valor pela Monotonicidade c/ Acelerador
tic;
[v4, idx4] = iteration_monotonicity_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(4) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx4);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com acelerador s�o %s. \n\n', teste);

%% Concavity
clc;

% M�todo 3: Itera��o da Fun��o Valor pela Concavidade
tic;
[v5, idx5] = iteration_concavity(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(5) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx5);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Concavidade s�o %s. \n\n', teste);

%% Concavity w/ Accelerator
clc;

% M�todo 3: Itera��o da Fun��o Valor pela Concavidade c/ Acelerador
tic; 
[v6, idx6] = iteration_concavity_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(6) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx6);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Concavidade com acelerador s�o %s. \n\n', teste);

%% Monotonicity & Concavity
clc;

% M�todo 4: Itera��o da Fun��o Valor pela Monotonicidade & Concavidadade
tic;
[v7, idx7] = iteration_monotonicity_and_concacity(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(7) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx7);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Concavidade com acelerador s�o %s. \n\n', teste);

%% Monotonicity & Concavity w/ Accelerator
clc;

% M�todo 4: Itera��o da Fun��o Valor pela Monotonicidade & Concavidade c/ Acelerador
tic;
[v8, idx8] = iteration_monotonicity_and_concacity_a(v, kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer(8) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx1,idx8);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Concavidade com acelerador s�o %s. \n\n', teste);

%% Tempos
Timer = timer';
Method = ["BF"; "BFwA"; "M"; "MwA"; "C"; "CwA"; "CnM"; "CnMwA"];
varTypes = ["string", "double"];
table(Method, Timer)

%% Multigrid
clc; clear v;

% Escolha do grid mais grosso ao grid mais fino
multigrid = [100, 500, 5000];
nk        = multigrid(end);
kgrid     = linspace(kmin, kmax, nk)';

% Benjamin Moll's guess (melhor que come�ar do zero)
v         = zeros(multigrid(1), nz);     
init_grid = linspace(kmin, kmax, multigrid(1))';    % inicio o guess com o grid mais grosso

for iz = 1:nz
    for ik = 1:multigrid(1)
        v(ik,iz) = u(zgrid(iz).*init_grid(ik).^alpha + (1-delta).*init_grid(ik) - init_grid(ik))./(1-beta); % esse ser� meu guess inicial
    end
end

% Iniciar o processo de multigrid

% Brute Force
tic;
[v91, idx91] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_bruteforce, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(1) = toc;

% Brute Force with Accelerator
tic;
[v92, idx92] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_bruteforce_a, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(2) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx92);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Brute Force com acelerador s�o %s. \n\n', teste);

% Monotonicity 
tic;
[v93, idx93] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_monotonicity, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(3) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx93);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade s�o %s. \n\n', teste);

% Monotonicity with Accelerator
tic;
[v94, idx94] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_monotonicity_a, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(4) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx94);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com Acelerador s�o %s. \n\n', teste);

% Concavity 
tic;
[v95, idx95] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_concavity, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(5) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx95);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com Acelerador s�o %s. \n\n', teste);

% Concavity with Accelerator
tic;
[v96, idx96] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_concavity_a, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(6) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx96);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com Acelerador s�o %s. \n\n', teste);

% Monotonicity and Concavity
tic;
[v97, idx97] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_monotonicity_and_concacity, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(7) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx97);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com Acelerador s�o %s. \n\n', teste);

% Monotonicity and Concavity with Accelerator
tic;
[v98, idx98] = iteration_multigrid(v, ... 
                                 multigrid, init_grid, ...
                                 @iteration_monotonicity_and_concacity_a, ...
                                 zgrid, P, alpha, beta, delta, mu, max_iter, tol);
timer_mg(8) = toc;

% As matrizes de �ndices s�o iguais?
igualdade = isequal(idx91,idx98);
if igualdade == 1; teste = 'iguais'; else; teste = 'diferentes'; end
fprintf('As matrizes de �ndices dos modelos de Brute Force e Monotonicidade com Acelerador s�o %s. \n\n', teste);

%%

% Decido usar o 97 por ser mais seguro e r�pido
% Extrair as fun��es pol�tica k' = g(k,z):
g97 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        g97(ik,iz) = kgrid(idx97(ik,iz));
    end
end

% Extrair as fun��es pol�tica c = c(k,z):
c97 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        c97(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha)+(1-delta)*kgrid(ik)-g97(ik,iz);
    end
end

%% Gr�fico - Fun��o Valor (2D & 3D)
figure(5)
plot_value_function(v97, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Capital) (2D & 3D)
figure(6)
plot_capital_policy_function(g97, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Consumo) (2D & 3D)
figure(7)
plot_consumption_policy_function(c97, kgrid, zgrid);

%% Euler Equations Error (EEE)
figure(8)
EEE2 = euler_equation_erros(idx97, c97, kgrid, zgrid, P, alpha, beta, delta, mu);

% amplitude
max(max(EEE2))
min(min(EEE2))

%% Grid end�geno 
clc;

% Definir novamente o grid para 500
nk       = 500; 
kmax     = 1.25*kss;
kmin     = 0.75*kss;

kgrid    = linspace(kmin, kmax, nk)';

% Iniciar o processo do grid end�geno
tic;
c10 = iteration_egm(kgrid, zgrid, P, alpha, beta, delta, mu, max_iter, tol);

% Fun��o pol�tica do capital (grid end�geno)
g10 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        g10(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - c10(ik,iz);
    end
end

% Fun��o pol�tica do capital (for�ando estar no grid)
g10_aprox = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        dif = abs(kgrid - g10(ik,iz));
        lb  = min(dif);
        g10_aprox(ik,iz) = kgrid(dif == lb);        
    end
end

% Fun��o pol�tica do consumo (for�ando estar no grid)
c10_aprox = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        c10_aprox(ik,iz) = zgrid(iz)*(kgrid(ik)^alpha) + (1-delta)*kgrid(ik) - g10_aprox(ik,iz);
    end
end
timer(11) = toc;

%% Gr�fico - Fun��o Pol�tica (Capital) (2D & 3D) - Grid end�geno
figure(9)
plot_capital_policy_function(g10, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Capital) (2D & 3D) - Grid ajustado
figure(10)
plot_capital_policy_function(g10_aprox, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Consumo) (2D & 3D) - Grid end�geno
figure(11)
plot_consumption_policy_function(c10, kgrid, zgrid);

%% Gr�fico - Fun��o Pol�tica (Consumo) (2D & 3D) - Grid ajustado
figure(12)
plot_consumption_policy_function(c10_aprox, kgrid, zgrid);

%% Recuperando a fun��o valor

% Recuperando a matriz de �ndices
idx10 = zeros(nk,nz);
for iz = 1:nz
    for ik = 1:nk
        idx10(ik,iz) = find(kgrid == g10_aprox(ik,iz));
    end
end

% Matrizes de aloca��o
% --------------------
iter  = 0;
error = 1;
v10   = zeros(nk,nz);            % Guess inicial;
Tv    = zeros(nk,nz);            % Matriz do operador de Bellman;

while (error > tol && iter <= max_iter)
    for iz = 1:nz        
        for ik = 1:nk
            
            % Pegando o capital ik e criando um vetor de consumo
            % -----------------------------------------------------------
            k = kgrid(ik);      % por algum motivo aumenta a performance do c�gido
           
            % n�o h� consumo �timo negativo
            u0 = (c10_aprox(ik,iz)^(1-mu)-1)/(1-mu);
            
            %%% Computar a esperan�a
            Ev = 0;
            for jz = 1:nz
                Ev = Ev + P(iz,jz)*v10(idx10(ik,iz),jz);
            end
            
            %%% Computar a fun�ao valor (matriz)
            Tv(ik,iz) = u0 + beta.*Ev;
        end
    end
    
    error = max(max(abs(Tv-v10)));
    iter  = iter + 1;                   % atualiza a  itera��o
    v10 = Tv;
    
    % Imprime a itera��o e o erro
    fprintf('Error %4i %6.2e \n',[iter, error]);
end

% Plot
figure(13)
plot_value_function(v10, kgrid, zgrid)

%% Euler Equation Erros (Grid End�geno)
figure(14)
EEE3 = euler_equation_erros_egm(c10, g10_aprox, kgrid, zgrid, P, alpha, beta, delta, mu);

max(max(EEE3))
min(min(EEE3))

%% Comparando os Erros de Euler
figure(15)
set(gcf,'DefaultLineLineWidth', 1.5);
plot(kgrid,EEE1(:,4),kgrid,EEE3(:,4))
     set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('$\log10$ $\mid$ Euler Equation Error $\mid$','Interpreter','latex');
title('\textbf{Equa\c{c}{\~a}o dos Erros de Euler (EEE) para $$z = 1$$}','interpreter','latex','fontsize',20);
ha = legend('Grid Search (Brute Force)', 'M{\''e}todo End{\''o}geno',...
            'Location','SouthEast');
set(ha,'Interpreter','latex');
axis([kmin kmax 1.10*min(min(EEE1)) 0.90*max(max(EEE1))])
