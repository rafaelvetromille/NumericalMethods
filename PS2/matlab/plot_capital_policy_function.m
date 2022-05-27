function plot_capital_policy_function(variable, kgrid, zgrid)

% Necessário para plotar o gráfico em 3D
k            = repelem(kgrid,1,7);
z            = repelem(zgrid,1,length(kgrid))';

% Gráfico 2 - Função política do capital
subplot(1,2,1);

% Gráfico (3D) - Função política: k' = g(k,z)
% -------------------------------------------
set(gcf,'DefaultLineLineWidth', 1.5);
surf(k, z, variable, 'FaceAlpha', 0.8, 'EdgeColor', 'none')
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('realiza\c{c}{\~o}es de produtividade, $z$','Interpreter','latex');
zlabel('fun\c{c}{\~a}o  pol{\''i}tica, $k^\prime = g(k,z)$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o pol{\''i}tica $$k^\prime = g(k,z)$$}','interpreter','latex','fontsize',20);
axis([min(kgrid) max(kgrid) min(zgrid) max(zgrid) min(min(variable)) max(max(variable))]) 

% Gráfico (2D) - Função política: k' = g(k,z)
% -------------------------------------------
subplot(1,2,2);
set(gcf,'DefaultLineLineWidth', 1.5);
plot(kgrid,variable)
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
set(get(gcf,'CurrentAxes'),'FontSize',14,'LineWidth', 1.5)
set(gca,'DefaultTextInterpreter', 'latex','TickLabelInterpreter','latex') 
xlabel('estoque de capital, $k$','Interpreter','latex');
ylabel('fun\c{c}{\~a}o pol{\''i}tica, $k^\prime = g(k,z)$','Interpreter','latex');
title('\textbf{Fun\c{c}{\~a}o pol{\''i}tica $k^\prime = g(k,z)$ para v{\''a}rios $$z$$''s}','interpreter','latex','fontsize',20);
ha = legend('$z=z_1$', '$z=z_2$', '$z=z_3$',...
            '$z=z_4$', '$z=z_5$', '$z=z_6$',...
            '$z=z_7$',...
            'Location','SouthEast');
set(ha,'Interpreter','latex');
axis([min(kgrid) max(kgrid) min(min(variable)) max(max(variable))]) 

end