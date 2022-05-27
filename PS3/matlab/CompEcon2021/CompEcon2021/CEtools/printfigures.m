%% PRINTFIGURES
%
% Save figures for overheads and manuscripts

%  Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function printfigures(nam)

return

h = findobj('type','figure');
for i=1:length(h)
  figure(i)
  title([])
  box off
  eval(['print -depsc ' ['..\Figures\' nam num2str(i)]])
end