function [Nv] = densityInValence(T,mhStar)
% This is calculating formula 3.3 in
% W. Shockley and W.T. Read, Jr. Phys. Rev. 87, 835 (1952)
% Statistics of the Recombinations of Holes and Electrons

%The integral is approximating fermi-Dirac statistics to Boltzmann
%statistics, the approximation may fail
% TODO: compare these calculated values with Fermi-Dirac statistics


%Boltzmann constant in Joules/Kelvin
kb = 1.380649e-23;

%BoltzmanEnergy in joules
kbT = kb.* T;

%Planck's constant in Joules-seconds
h = 6.62607015e-34;
hbar = h/(2*pi);

%Electron rest mass
me0 = 9.10938370e-31;


%mhStar  =   effective mass of the hole at the highest point of the valence
%band
%Putting it in units of cm^-3
Nv = 1/4.* (2.*mhStar*me0.*kbT/(pi*hbar^2)).^(3/2)./(100)^3;

