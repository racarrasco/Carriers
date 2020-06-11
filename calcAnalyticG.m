function [radiativeGenCoefficient] = calcAnalyticG(x, meStar, mhStar, ni, Eg, einf)
% Calculate the equilibrium rate of radiative recombination found from
% references M. A. Kinch, M. J. Brau, And A. Simmons, JAP vol 44, p. 1649
% (1973) 
% https://doi.org/10.1063/1.1662426
%
% This could also be found in R. N. Hall  Proc. Inst. of Electr. Eng. B vol
% 106 p 923 (1959).Recombination Processes in Semiconductors.
% 10.1049/pi-b-2.1959.0171



%Make the temperature as a row because the band gaps are a row
if(size(x,1) > 1)
    T = x';
else
    T = x; % the temperature is already a row
end



%Free electron mass in kg
meFree = 9.10938356e-31;



%{ 
Follows the formula from eq
radiativeGenCoefficient = ni.^2 .* 5.8e-13 .* einf^(1/2) .*...
    (meFree/(meStar.*meFree + mhStar.*meFree))^(3/2).* ...
    (1 + 1/meStar + 1/mhStar) .* (300./T).^(3/2).*Eg.^2;
%}





% boltzmann constant in joules/Kelvin
boltzmann = 1.380649e-23;

% Electron charge in Coulomb
e = 1.602176634e-19;  

% Boltzmann constant in eV/Kelvin
keV = boltzmann/e;


% Thermal energy in Joules
kT = keV.*T;

% Calculate according to V. C. Lopes, A. J. Syllaios and M. C. Chen,
% Semicond. Sci. Technol. vol. 8, p. 824 (1993). Minority carrier lifetime
% in mercury cadmium telluride. 
% doi:  [10.1088/0268-1242/8/6S/005]

radiativeGenCoefficient = ni.^2 .* 5.8e-13 .* einf^(1/2) .*...
    (meFree/(meStar.*meFree + mhStar.*meFree))^(3/2).* ...
    (1 + 1/meStar + 1/mhStar) .* (300./T).^(3/2).* ...
    (Eg.^2 + 3.*kT.*Eg + 3.75.*(kT).^2);


