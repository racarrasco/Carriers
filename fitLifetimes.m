function[fitObject, gof, output] = fitLifetimes(... 
 xin,y,  type, meStar, mhStar, eg, eg300, valenceEdge, f1f2, conductionEdge,...
 einf, nRefractive, kExtinction, photonEnergyum, t, n, dopingDensity, defectLevel, defectDensity)


% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% SEE THE calculateLifetimes.m file for meaning to the input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xin =  Temperature array that goes up to the resolution
% y   = results of the extracted TRPL data
% a   = type of doping, n-type, p-type or intrinsic
% b   = meStar electron effective mass in units of electron rest masses
% c   = mhStar, heavy hole effective mass in units of electron rest masses
% d   = eg, band gap as a function temperature
% e   = eg300 the band gap at 300 K
% f   = Ev valence band edge
% g   = f1f2 overlap integral
% h   = conduction edge minimum
% k  = static refractive index einf (for auger)
% l  = refractive index as a function of photon energy
% m   = extinction coefficient as a function of photon energy
% n   = photon energy in micron
% o   = active layer thickness for photon recycling
% p   = refractive index of the cap layer
% q   = doping density
% r   = defect level
% s   = defect density

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, not in the order in which they appear.

guess = [dopingDensity, defectLevel, defectDensity];

% Only fit 3 parameters which are o, p, and q
ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s)',...
    'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o'; 'p']);

%The three parameters will be the doping density, 
[fitObject, gof, output] = fit(xin',y',ft,'problem',{type, meStar, mhStar, eg, eg300, ...
    valenceEdge, f1f2, conductionEdge, einf,nRefractive,kExtinction,photonEnergyum, t,n}, ...
    'Startpoint', guess);