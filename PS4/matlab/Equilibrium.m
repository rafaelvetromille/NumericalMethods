clear

%%

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

if j == 1

termo = (z_1 - rho*grid_tauchen(i) + delta_z/2)/sqrt(sigma_sq);

p_i_j = normcdf(termo);

end

%%%%%%%%%%%%% MAIOR %%%%%%%%%%%%;

if j == N

termo = (z_N - rho*grid_tauchen(i) - delta_z/2)/sqrt(sigma_sq);

p_i_j = 1 - normcdf(termo);

end

P(i,j) = p_i_j;

end

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% CRIAR O GRID DOS ATIVOS %%%%%%%%%%%%

grid_z = grid_tauchen;

%%%%%%%%%%%%%%%%%%%%%% LIMITANDO a e b

tic

excessos_de_demanda = 0;

for i = 1:8

excessos_de_demanda(i) = excess_demand(i/200);

i

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USANDO MONOTONICIDADE E CONCAVIDADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% VAMOS RESOLVER PELO MÉTODO DA BISSETRIZ

b = 0.9999*(1/beta - 1);
a = 0.04;

%%%% J�? CHECAMOS QUE COM b A DEMANDA É POSITIVA E COM a É NEGATIVA


erro_bissetriz = 1;

while(erro_bissetriz > tol)

c = (a + b)/2;

demanda_a = excess_demand(a);
demanda_b = excess_demand(b);
demanda_c = excess_demand(c);

if demanda_a*demanda_c < 0

b = c;

end

if demanda_b*demanda_c < 0

a = c;

end

erro_bissetriz = abs(demanda_c - 0)

end

time_elapsed = toc