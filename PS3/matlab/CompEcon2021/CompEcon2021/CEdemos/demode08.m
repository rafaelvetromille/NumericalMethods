%% DEMODE08 Lorentz Strange Attractor

%  Solve
%   x1' = -a*x1 +a*x2
%   x2' = -x2 - x1*x3
%   x3' = -b*x3 + x1*x2 - b*c;
%  s.t.
% 	x1(0) = -8
% 	x2(0) =  8
% 	x3(0) = 27

% Preliminary tasks
deminit(mfilename)


%% FORMULATION

% Velocity Function
a = 10;
b = 8/3;
c = 28;
f = @(x) [-a*x(1,:)+a*x(2,:);-x(2,:)-x(1,:).*x(3,:);-b*x(3,:) + x(1,:).*x(2,:)-b*c];

% Initial State
xinit = [-8;8;27];
               
% Time Horizons
T = 12;


%% SOLVE ODE USING RUNGE-KUTTA METHOD (ODERK4)

% Solve ODE
N = 750;   % number of time nodes
[x,t] = oderk4(f,xinit,T,N);

% Plot Solution
figure
plot(t,x)
title('Lorentz Strange Attractor')
xlabel('Time')
ylabel('State')
legend('$x_1$','$x_2$','$x_3$')


%% PHASE DIAGRAM

% Layout
x1lim = [-20 20];       % x1 plotting limits
x2lim = [-30 30];       % x2 plotting limits
odephase(x1lim,x2lim, ...
  'Lorentz Strange Attractor Phase Diagram','$x_1$','$x_3$')

% Plot State Path
hold on
for j=2:N
   plot(x(j-1:j,1),x(j-1:j,3),'r','LineWidth',2)
   getframe; 
end
plot(x(:,1),x(:,3),'r','LineWidth',2)


%% SAVE FIGURES
printfigures(mfilename)