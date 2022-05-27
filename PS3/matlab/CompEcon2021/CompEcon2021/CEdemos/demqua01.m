%% DEMQUA01 Equidistributed Sequences on Unit Square

% Preliminary tasks
deminit(mfilename)

% Number of integration nodes
n = 2500;

% Neiderreiter Sequence
figure
x = qnwequi(n,[0 0],[1 1],'N');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
axis square
xticks([0 1])
yticks([0 1])
title('Neiderreiter Sequence')
xlabel('$x_1$')
ylabel('$x_2$')

% Weyl Sequence
figure
x = qnwequi(n,[0 0],[1 1],'W');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
axis square
xticks([0 1])
yticks([0 1])
title('Weyl Sequence')
xlabel('$x_1$')
ylabel('$x_2$')

% Pseudo-Random Sequence
figure
x = qnwequi(n,[0 0],[1 1],'R');
plot(x(:,1),x(:,2),'.','MarkerSize',7);
axis square
xticks([0 1])
yticks([0 1])
title('Pseudo-Random Sequence')
xlabel('$x_1$')
ylabel('$x_2$')


%% SAVE FIGURES
printfigures(mfilename)