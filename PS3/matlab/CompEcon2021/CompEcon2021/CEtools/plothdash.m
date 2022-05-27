%% PLOTHDASH
%
% Plots horizontal dashed line from y-axis to point (x,y) on current figure
%
%  Usage
%    plothdash(x,y,color,width)
%  Input
%    x      : x-coordinate of point
%    y      : y coordinate of point
%    v      : color of line (optional, default is black)
%    width  : width of line (optional, default is 2)
%
%  Copyright(c) 2016-2021
%    Mario J. Miranda - miranda.4@osu.edu

function plothdash(x,y,color,width)
xl = xlim;
if isempty(x), x=xl(2); end 
if nargin<3||isempty(color), color='k'; end
if nargin<4||isempty(width), width=2; end
plot([xl(1) x],[y y],':','LineWidth',width,'color',color)