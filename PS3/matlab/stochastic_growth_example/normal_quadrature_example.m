%%%%%  normal quadrature example

clearvars

%%%%% x is N(mu,sigma^2)

n = 3;        %% nodes

sigma = 1;    %% std dev of log(x)
mu    = 0;    %% mean of log(x)

[x,w] = qnwnorm(n,mu,sigma^2);


figure(1)
plot(x,w)
xlabel('quadrature node')
ylabel('quadrature weight')

fprintf('\n')
fprintf('quadrature nodes and weights \n')
fprintf('\n')
fprintf('  nodes       = %9.3f\n', x)
fprintf('  weights     = %9.3f\n' ,w)
fprintf('\n')



