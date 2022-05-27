%% Contents
%
% Demonstration programs to accompany
%   Lecture Notes in Computational Economic Dynamics
%   Professor Mario J. Miranda  
%   The Ohio State University
%
% Requires installation of
%   Matlab version 2019b or later
%   CompEcon2021 toolbox available for download at
%   https://aede.osu.edu/our-people/mario-javier-miranda
%   
% Based on book
%   Applied Computational Economics and Finance
%   Mario J. Miranda & Paul L. Fackler
%   2002, MIT Press, Cambridge MA
%
% Programs are regularly revised
%   This version 4/11/2021
%
% Copyright(c) 1997-2021
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

close all
clear all

% Introduction 
demintro01 % Inverse Demand Problem
demintro02 % Rational Expectations Agricultural Market Model
 
% Mathematics Review 
demmath01  % Taylor Approximations
demmath02  % Function Inner Products, Norms & Metrics
demmath03  % Discrete and Continuous Distributions
demmath04  % Standard Copulas
demmath05  % Inverse Function, Implicit Function and Mean Value Theorems
demmath06  % Operations with Markov Chains
 
% Linear Equations 
demlin01   % Linear Equations 
demlin02   % Ill-Conditioning of Vandermonde Matrices
demlin03   % Sparse Linear Equations
 
% Nonlinear Equations 
demslv01   % Root of f(x)=exp(-x)-1
demslv02   % Root of Rosencrantz Function
demslv03   % Fixedpoint of f(x) = x^0.5
demslv04   % Fixedpoint of g(x1,x2)= [x1^2+x2^3;x1*x2-0.5]
demslv05   % Newton and Broyden Paths
demslv06   % Nonlinear Equation Methods
demslv07   % Linear Complementarity Problem Methods
demslv08   % Nonlinear Complementarity Problem
demslv09   % Hard Nonlinear Complementarity Problem with Billup's Function
demslv10   % Linear Complementarity Problem
demslv11   % Rootfinding Reformulations of Nonlinear Complementarity Problem
demslv12   % Convergence Rates for Nonlinear Equation Methods
demslv13   % Spatial Equilibrium Model
 
% Finite-Dimensional Optimization
demopt01   % Maximization via Golden Search
demopt02   % Changes in Nelder-Mead Simplex
demopt03   % Nelder-Mead Simplex Method
demopt04   % Maximization of Rosencrantz Function by Various Methods
demopt05   % Optimization with qnewton
demopt06   % KKT Conditions for Constrained Optimization Problems
demopt07   % Bound-Constrained Optimization via Sequential LCP
demopt08   % Constrained Optimization with nlpsolve
demopt09   % Linear Programming lpsolve
 
% Quadrature
demqua01   % Equidistributed Sequences on Unit Square
demqua02   % Expectation of Function of Random Normal Vector
demqua03   % Area Under a Curve, Various Methods
demqua04   % Area under Normal PDF Using Simpson's Rule
demqua05   % Expected Utility Insurance Model
demqua06   % Area Under a Curve
demqua07   % Monte Carlo Simulation of Time Series

% Numerical Differentiation
demdif01   % Finite Difference Hessian Evaluation Structure
demdif02   % Error in Finite Difference Derivatives
demdif03   % Demonstrates fdjac and checkjac
demdif04   % Demonstrates fdhess
demdif05   % Finite-Difference Jacobians and Hessians
 
% Function Approximation
demapp01   % Approximating Functions on R
demapp02   % Approximating Functions on R^2
demapp03   % Standard Basis Functions and Nodes
demapp04   % Uniform- and Chebychev-Node Polynomial Approximation of Runge's Function
demapp05   % Chebychev Polynomial and Spline Approximation of Various Functions
demapp06   % Cournot Oligopoly Model
demapp07   % Compute Implicit Function via Collocation
demapp08   % Linear Spline Approximation
demapp09   % Monopolist's Effective Supply
 
% Discrete Time Discrete State Dynamic Programming
demddp01   % Mine Management Model
demddp02   % Asset Replacement Model
demddp03   % Binomial American Put Option Model
demddp04   % Water Management Model
demddp05   % Bioeconomic Model
demddp06   % Renewable Resource Model
demddp07   % Job Search Model
demddp08   % Deterministic Cow Replacement Model
demddp09   % Stochastic Cow Replacement Model
demddp10   % Stochastic Optimal Growth Model

% Discrete Time Continuous State Dynamic Programming
demdp00    % Miscellaneous Lecture Note Figures
demdp01    % Timber Harvesting Model
demdp02    % Asset Replacement Model
demdp03    % Industry Entry-Exit Model
demdp04    % Job Search Model
demdp05    % American Put Option Pricing Model
demdp06    % Ramsey Stochastic Optimal Economic Growth Model
demdp07    % Stochastic Economic Growth Model
demdp08    % Public Renewable Management Model
demdp09    % Private Non-Renewable Resource Management Model
demdp10    % Water Resource Management Model
demdp11    % Monetary Policy Model
demdp12    % Production Management Model
demdp13    % Inventory Management Model
demdp14    % Livestock Feeding Model
demdp15    % Lifecycle Consumption-Savings Model
demdp16    % Savings and Insurance Model
demdp17    % Credit With Technology Adoption Model
demdp18    % Linear-Quadratic Model
demdp19    % Credit With Strategic Default Model
demdp20    % Saving with Transactions Costs Model

% Discrete-Time Rational Expectations Models
demrem01   % Asset Pricing Model
demrem02   % Commodity Storage Model
demrem03   % Government Price Support Model
 
% Discrete-Time Dynamic Game Models
demgame01  % Production Capacity Game Model
demgame02  % Income Redistribution Game Model
demgame03  % Marketing Board Game Model

% Ordinary Differential Equations
demode01   % Stability of Linear ODEs
demode02   % Nonlinear ODE Example
demode03   % Linear ODE Example
demode04   % Non-IVP Linear ODE Example
demode05   % Commodity Storage Model
demode06   % Predator-Prey Model
demode07   % Commercial Fisheries Model
demode08   % Lorentz Strange Attractor
demode09   % Tobin's Q Model
demode10   % Regional Migration Model

% Continuous Time Deterministic Optimal Control
demdoc01   % Optimal Consumption-Investment Model
demdoc02   % Optimal Economic Growth Model
demdoc03   % Nonrenewable Resource Model
demdoc04   % Renewable Resource Model
demdoc05   % Production Adjustment Model

% Continuous Time Stochastic Optimal Control
demsoc00   % Ito Processes
demsoc01   % Consumption-Investment Model
demsoc02   % Portfolio Model
demsoc03   % Optimal Economic Growth Model
demsoc04   % Renewable Resource Model
demsoc05   % Optimal Fish Harvest Model
demsoc06   % Production-Adjustment Model
demsoc07   % Nonrenewable Resource Model

% Regime Switching Models
demrs01    % Asset Abandonment Model
demrs02    % Fish Harvest Model
demrs03    % Dixit Entry/Exit Model
 
% Impulse Control Models
demic01    % Asset Replacement Model
demic02    % Timber Harvesting Model
demic03    % Storage Management Model
demic04    % Capacity Choice Model
demic05    % Cash Management Model
demic06    % Fish Harvest Model