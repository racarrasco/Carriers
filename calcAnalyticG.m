function [radiativeGenCoefficient] = calcAnalyticG(T, meStar, mhStar, ni, Eg)
% Calculate the equilibrium rate of radiative recombination found from
% references M. A. Kinch, M. J. Brau, And A. Simmons, JAP vol 44, p. 1649
% (1973) 
% https://doi.org/10.1063/1.1662426
%
% This could also be found in R. N. Hall  Proc. Inst. of Electr. Eng. B vol
% 106 p 923 (1959).Recombination Processes in Semiconductors.
% 10.1049/pi-b-2.1959.0171

%{
% boltzmann constant in joules/Kelvin
boltzmann = 1.380649e-23;

% Electron charge in Coulomb
e = 1.602176634e-19;  

% Boltzmann constant in eV/Kelvin
keV = boltzmann/e;

% Thermal energy in Joules
kT = keV.*T;
%}


radiativeGenCoefficient = ni.^2 .* 5.8e-13 .* einf^(1/2) .*...
    (m0/(meStar.*m0 + mhStar.*m0))^(3/2).* ...
    (1 + 1/meStar + 1/mhStar) .* (300./T).^(3/2).*Eg.^2; %...
   % (Eg.^2 + 3.*kT.*Eg + 3.75.*(kT).^2);
