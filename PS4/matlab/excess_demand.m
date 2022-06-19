function F = excess_demand(r)



%%%%%%%%%% DEFININDO PARÂMETROS %%%%%%%

m = 3;
sigma_sq = 0.01^2;
rho = 0.9;

%%% TAMANHO DO GRID DOS CHOQUES %%%%%%%

N = 9;

beta = 0.96;
mu = 1.0001;


sigma_sq_z = sigma_sq/(1-rho^2);

%%%%%%% TAMANHO DO GRID DOS ATIVOS

n_a = 500;

%%%%%%% TOLERANCIA

tol = 10^(-6);

tol_distrb = 10^(-5);

%%%%%%%% TAUCHEN %%%%%%%%%%%



%%%%%%% CRIANDO O grid_tauchen %%%%%%%%%%%%%;

z_N = m*sqrt(sigma_sq_z);

z_1 = -z_N;

delta_z = (z_N - z_1)/(N-1);

grid_tauchen = z_1:delta_z:z_N;


%%%%%%%%%%%%%%%%% PREENCHENDO A MATRIZ P PARA TAUCHEN %%%%%%%%%%%%;

for i = 1:N 

for j = 1:N


termo_1 = (grid_tauchen(j) + delta_z/2 - rho*grid_tauchen(i))/sqrt(sigma_sq);

termo_2 = (grid_tauchen(j) - delta_z/2 - rho*grid_tauchen(i))/sqrt(sigma_sq);

p_i_j = normcdf(termo_1) - normcdf(termo_2);

%%%%%%%%%%%%% CORRIGINDO CORNERS %%%%%%%%;

 %%%%%%%%%%%%% MENOR %%%%%%%%%%;

if j == 1;

termo = (z_1 - rho*grid_tauchen(i) + delta_z/2)/sqrt(sigma_sq);

p_i_j = normcdf(termo);

end

%%%%%%%%%%%%% MAIOR %%%%%%%%%%%%;

if j == N;

termo = (z_N - rho*grid_tauchen(i) - delta_z/2)/sqrt(sigma_sq);

p_i_j = 1 - normcdf(termo);

end

P(i,j) = p_i_j;

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% CRIAR O GRID DOS ATIVOS %%%%%%%%%%%%


grid_z = grid_tauchen;


%%%%%%%%%%% PREENCHER A V INICIAL %%%%%%%%%%%%%%


%%%%%%% ORDEM DA FUNÇÃO VALOR: LINHA É a, COLUNA É Z %%%%%%%%%%

V_0 = zeros(n_a, N);

r_max = 1/beta - 1;

limite_natural = -exp(grid_z(1))/r_max;

grid_a = linspace(limite_natural, -limite_natural, n_a);


V_new = V_0;



erro = 1;

while(erro > tol)

V_old = V_new;

for i = 1:N

probabilidades = P(i,:);

for j = 1:n_a

estoque = exp(grid_z(i)) + (1+r)*grid_a(j);

%%%%%%% START = 1 CASO j = 1

if j == 1;

start = 1;

end

%%%%%% CALCULAR OS VALORES EM TODOS OS a' POSS�?VEIS %%%%%%%%%%%

valor = zeros(1,n_a);

%%% PRIMEIRO, COMPUTAR PARA t = 1

t = start;


c_t = estoque - grid_a(t);

%%% CALCULAR O VALOR DA UTILIDADE, CASO C_T < 0, U = -INFINITO

if c_t >= 0;

u = c_t^(1-mu) - 1;
u = u/(1-mu);

elseif c_t < 0;

u = -Inf;

%%% ENCERRAR IF

end 

%%%%%%%% CALCULAR A ESPERANÇA

esperanca = 0;

for w = 1:N

esperanca = esperanca + probabilidades(w)*V_old(t,w);

%%% ENCERRAR O FOR DE w

end


valor(t) = u + beta*esperanca;



%%%%%%%%%%%%%%%%% FAZER O FOR A PARTIR DE START + 1

for t = (start+1):n_a

c_t = estoque - grid_a(t);

%%% CALCULAR O VALOR DA UTILIDADE, CASO C_T < 0, U = -INFINITO

if c_t >= 0;

u = c_t^(1-mu) - 1;
u = u/(1-mu);

elseif c_t < 0;

u = -Inf;

%%% ENCERRAR IF

end 

%%%%%%%% CALCULAR A ESPERANÇA

esperanca = 0;

for w = 1:N

esperanca = esperanca + probabilidades(w)*V_old(t,w);

%%% ENCERRAR O FOR DE w

end


valor(t) = u + beta*esperanca;

%%%%%%%%%%%%%%%%%%%% AQUI EST�? A DIFERENÇA: USAMOS A CONCAVIDADE

if valor(t) < valor(t-1);

break

end

%%% ENCERRAR O FOR DE t

end

%%%%%%%%%%%%%%%%%%%%%% VEJA QUE NÃO H�? MAX

%%%%% MUDANÇA CASO NUNCA TENHA ACIONADO O BREAK (OU SEJA, O M�?XIMO É NO ÚLTIMO t MESMO)

if valor(t) < valor(t-1);

V_new(j,i) = valor(t-1);

elseif valor(t) >= valor(t-1);

V_new(j,i) = valor(t);

end

%%%%% MUDANÇA CASO NUNCA TENHA ACIONADO O BREAK (OU SEJA, O M�?XIMO É NO ÚLTIMO t MESMO)

if valor(t) < valor(t-1);

indice_a = t-1;

elseif valor(t) >= valor(t-1);

indice_a = t;

end

%%%%%%%%%%%%%%%%%%%% AQUI EST�? OUTRA DIFERENÇA: NOSSO START É A PARTIR DA DECISÃO DO INDIV�?DUO ANTERIOR (MONOTONICIDADE)

start = indice_a;

g(j,i) = grid_a(indice_a);


%%% ENCERRAR FOR DE j

end

%%% ENCERRAR FOR DE i

end

%%%%%%%%%%%

%%%% CALCULAR ERRO

dif = abs(V_new - V_old);

erro = max(max(dif));

%%%%%%% ENCERRAR WHILE


end





%%%%%%%%%%%%%% CALCULAR A DISTRIBUIÇÃO INVARIANTE

%%% CHUTE INICIAL pi_0

pi_0 = ones(n_a, N)/(N*n_a);

pi_new = pi_0;

erro = 1;

while(erro > tol_distrb)

pi_old = pi_new;

%%%%%% ITERAR

for i = 1:N
for j = 1:n_a

select = (g == grid_a(j));

aux = pi_old.*select;


pi_new(j,i) = sum(aux)*P(:,i);

end
end

%%%% CALCULAR ERRO

dif = abs(pi_new - pi_old);

erro = max(max(dif));

%%%%%%% ENCERRAR WHILE


end

%%% CALCULAR DEMANDA POR CAPITAL

demanda = 0;

for i = 1:N
for j = 1:n_a

demanda = demanda + g(j,i)*pi_new(j,i);

end
end

F = demanda;

end
