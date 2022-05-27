%% PLOTVDASH
%
% Plots vertical dashed line from x-axis to point (x,y) on current figure
%
%  Usage
%    plotvdash(x,y,color,width)
%  Input
%    x      : x-coordinate of point
%    y      : y coordinate of point
%    v      : color of line (optional, default is black)
%    width  : width of line (optional, default is 2)
%
%  Copyright(c) 2016-2021
%    Mario J. Miranda - miranda.4@osu.edu

function plotvdash(x,y,color,width)
yl = ylim;
if isempty(y), y=yl(2); end 
if nargin<3||isempty(color), color='k'; end
if nargin<4||isempty(width), width=2; end
plot([x x],[yl(1) y],':','LineWidth',width,'color',color)