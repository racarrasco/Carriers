function[fitObject, gof, output] = ...
 fitLifetimes(xin,y,  type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
 f1f2, einf, Nc, Nv, ni,G, phi,  dopingDensity, defectLevel, defectDensity)


% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% SEE THE calculateLifetimes.m file for meaning to the input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, not in the order in which they appear.

guess = [dopingDensity, defectLevel, defectDensity];

% Only fit 3 parameters which are o, p, and q
ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
    'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o']);





%{
fo = fitoptions('Lower', [0.000001, -1, 0],...
                'Upper', [10000,0, inf],...
                'StartPoint', guess);
%}
    
%The three parameters will be the doping density, defect level and defect
%density
[fitObject, gof, output] = fit(xin',y',ft,...
    'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge, f1f2,...
    einf, Nc, Nv, ni, G, phi},...
    'Lower', [0.000001, -1, 0],...
    'Upper', [100000, 0, inf], ...
    'Startpoint', guess);