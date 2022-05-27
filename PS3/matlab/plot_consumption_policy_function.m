function plot_consumption_policy_function(variable, kgrid, zgrid)

% Necessário para plotar o gráfico em 3D
k            = repelem(kgrid,1,7);
z            = repelem(zgrid,1,length(kgrid))';
nz           = length(zgrid);


% Gráfico (3D) - Função política: k' = g(k,z)
hold on
subplot(1,2,1);
set(gcf,'DefaultLineLineWidth', 1.5);
surf(k, z, variable, 'FaceAlpha', 0.8, 'EdgeColor', 'none')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('realiza\c{c}{\~o}es de produtividade, $z$','Interpreter','latex');
zlabel('fun\c{c}{\~a}o  pol{\''i}tica de consumo, $c = c(k,z)$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o pol{\''i}tica de consumo $$c = c(k,z)$$}','interpreter','latex','fontsize',20);
axis([min(kgrid) max(kgrid) min(zgrid) max(zgrid) min(min(variable)) max(max(variable))]) 
box on
hold off

% Gráfico (2D) - Função política: k' = g(k,z)
subplot(1,2,2);
set(gcf,'DefaultLineLineWidth', 1.5);
hold on
for i = 1:nz
    plot(kgrid, variable(:,i), 'DisplayName', sprintf('$z_{%4i} = %6.2f$', [i, zgrid(i)]));   
end
hl = legend('show', 'Location', 'southeast');
set(hl, 'Interpreter','latex')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex', 'TickLabelInterpreter','latex') 
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('fun\c{c}{\~a}o pol{\''i}tica de consumo, $c = c(k,z)$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o pol{\''i}tica de consumo $c = c(k,z)$ para v{\''a}rios $$z$$''s}','interpreter','latex','fontsize',20);
axis([min(kgrid) max(kgrid) min(min(variable)) max(max(variable))]) 
box on 
hold off

end