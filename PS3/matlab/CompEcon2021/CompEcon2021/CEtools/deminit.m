function deminit(filename)
disp(' ')
eval(['help ' filename])
clear all
close all
set(0,'defaultLineLineWidth',3) 
set(0,'defaultTextFontSize',12) 
set(0,'defaultAxesFontSize',14) 
set(0,'defaultLegendInterpreter','latex') 
set(0,'defaultTextInterpreter','latex')
set(0,'defaultAxesBox','off')
set(0,'defaultLegendBox','off')
set(0,'defaultLegendLocation','best')
set(0,'defaultTextHorizontalAlignment','left')
set(0,'defaultTextVerticalAlignment','bottom')
format short