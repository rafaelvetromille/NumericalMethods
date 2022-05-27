%%%%% lognormal quadrature example

clearvars

%%%%% log(x) is N(mu,sigma^2)

n = 11;       %% nodes

sigma = 0.2;  %% std dev of log(x)
mu    = 0;    %% mean of log(x)

[x,w] = qnwlogn(n,mu,sigma^2);

EX_quad   = sum(x.*w);
VarX_quad = sum((x.^2).*w)-sum(x.*w)^2;

%%%%% true moments

EX     = exp(mu+0.5*sigma^2);
VarX   = (exp(sigma^2)-1)*exp(2*mu+sigma^2);
xfine    = (0.01:0.01:5.00)';
density  = normpdf(log(xfine),mu,sigma)./xfine;

%%%%% sample moments

T = 1000;

logx = mu+sigma*randn(T,1);

EX_sim   = mean(exp(logx));
VarX_sim = var(exp(logx));

figure(1)
plot(x,w)
xlabel('quadrature node')
ylabel('quadrature weight')

fprintf('\n')
fprintf('mean \n')
fprintf('\n')
fprintf('  exact       = %9.3f\n', EX)
fprintf('  quadrature  = %9.3f\n' ,EX_quad)
fprintf('  simulated   = %9.3f\n' ,EX_sim)
fprintf('\n')

fprintf('\n')
fprintf('variance \n')
fprintf('\n')
fprintf('  exact       = %9.3f\n', VarX)
fprintf('  quadrature  = %9.3f\n' ,VarX_quad)
fprintf('  simulated   = %9.3f\n' ,VarX_sim)
fprintf('\n')


