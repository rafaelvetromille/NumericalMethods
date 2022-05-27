%% DEMODE02 Nonlinear Initial Value Problem

%  Solve
%    x1' = x1^2 - 2x2 - a
%    x2' = b - x1 - x2

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Phase Diagram Title and Axis Labels
figtitle = 'Nonlinear ODE Phase Diagram';
x1label  = '$x_1$';
x2label  = '$x_2$';

% Velocity Function
a =  5;
b = -1;
f = @(x) [x(1,:).^2-2*x(2,:)-a; b-x(1,:)-x(2,:)];

% Initial States
xinit = [-1 -1 -1; -3 -6 -4.5578];
       
% Time Horizons
T = [8 2 2];

       
%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Solve for Different Initial Values
N = 1000;   % number of time nodes
x = zeros(N,2,3);
for i=1:3
  [x(:,:,i),t] = oderk4(f,xinit(:,i),T(i),N);
end

% Plot Solutions in Time Domain
for i=1:3
  figure
  plot(t,x(:,:,i))
  xtickformat('%.1f')
  legend(x1label,x2label)
end


%% STEADY STATES

% Compute Steady State A
A1 = (-2-sqrt(4+4*(2*b+a)))/2;
A  = [A1;b-A1];
disp('Steady State A')
disp(A)
disp('Eigenvalues')
disp(eig(fdjac(f,A)))

% Compute Steady State B
B1 = (-2+sqrt(4+4*(2*b+a)))/2;
B  = [B1;b-B1];
disp('Steady State B')
disp(B)
disp('Eigenvalues')
disp(eig(fdjac(f,B)))


%% PHASE DIAGRAM

% Ploting Limits
x1lim = [-5 4];
x2lim = [-6 4];

% Compute Separatrix
T = 10;
xspx = odespx(f,B,T);

% Compute Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = [x1 (x1.^2-a)/2];
x2null = [x1 b-x1];

% Plot Phase Diagram
odephase(x1lim,x2lim,figtitle,x1label,x2label,f,x,xspx,[A B],x1null,x2null)

% Label State Paths
for i=1:3
  text(x(1,1,i)-0.1,x(1,2,i)+0.1,['\bf' num2str(i)])
end

% Label Steady States
text(A(1)+0.2,A(2)+0.1,'A')
text(B(1)-0.1,B(2)+0.2,'B')


%% SAVE FIGURES
printfigures(mfilename)