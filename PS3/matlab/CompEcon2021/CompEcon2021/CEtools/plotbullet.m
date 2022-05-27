%% PLOTBULLET
%
%  Inserts centered bullet in figure using simplified format.
%
%  Usage
%    bullet(x,y,size,color)
%  Input
%    x     : x coordinate
%    y     : y coordinate
%    size  : size of bullet (optional, default is 18)
%    color : color of bullet (optiional, default is black)
%
%  Copyright(c) 2016-2021
%    Mario J. Miranda - miranda.4@osu.edu

function plotbullet(x,y,size,color)
if nargin<3||isempty(size), size=18; end
if nargin<4||isempty(color), color='k'; end
if isempty(x)
  xl = xlim;
  x = xl(1);
end
if isempty(y)
  yl = ylim;
  y = yl(1);
end
text(x,y,'$\bullet$','color',color,'VerticalAlignment','middle','HorizontalAlignment','center','Fontsize',size,'Interpreter','Latex')