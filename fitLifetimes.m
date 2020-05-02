function[fitObject, gof, output] = ...
 fitLifetimes(xin,yin, fitAugerOverlap, type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
  einf, Nc, Nv, ni,G, phi, f1f2, dopingDensity, defectLevel, defectDensity)


% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% SEE THE calculateLifetimes.m file for meaning to the input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%yin   =  this will be a 1 x n x m matrix, where n is the length of the
%temperature array in which the measurements were performed, and m is the
%number of excitations performed. However, in order to perform fitting,the
%input needs to by and n x m matrix

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, not in the order in which they appear.

guess = [dopingDensity, defectLevel, defectDensity];
if(fitAugerOverlap) 
    ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
    'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n']);
   
    [fitObject, gof, output] = fit(xin',yin',ft,...
    'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
    einf, Nc, Nv, ni, G, phi},...
    'Lower', [0.01 , 0.00001, -1, 0],...
    'Upper', [0.3 , 100000, 0, inf], ...
    'Startpoint', guess);  
else
    % Only fit 3 parameters which are o, p, and q
ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
    'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o']);

[fitObject, gof, output] = fit(xin',yin',ft,...
    'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
    einf, Nc, Nv, ni, G, phi, f1f2},...
    'Lower', [0.000001, -1, 0],...
    'Upper', [100000, 0, inf], ...
    'Startpoint', guess);

end










%{
fo = fitoptions('Lower', [0.000001, -1, 0],...
                'Upper', [10000,0, inf],...
                'StartPoint', guess);
%}
    
%The three parameters will be the doping density, defect level and defect
%density

% reshape the input to a two-dimensional matrix
%y = reshape(yin,[size(xin,2), size(yin,3)]);
