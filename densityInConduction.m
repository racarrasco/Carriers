function [Nc] = densityInConduction(T,meStar)
% This is calculating formula 3.3 in
% W. Shockley and W.T. Read, Jr. Phys. Rev. 87, 835 (1952)
% Statistics of the Recombinations of Holes and Electrons

%The integral is approximating fermi-Dirac statistics to Boltzmann
%statistics, the approximation will fail for high doping and metals
%(degenerate electron gas)
%TODO: compare these calculated values with Fermi-Dirac statistics
%NOTE: we are working in SI units here, so energy is in Joules, not eV


%Boltzmann constant in Joules/Kelvin
kb = 1.380649e-23;

%BoltzmanEnergy in Joules
kbT = kb.* T;

%Planck's constant in Joules-seconds
h = 6.62607015e-34;
hbar = h/(2*pi);

%electron rest mass
me0 = 9.10938370e-31;


%meStar  =   effective mass of the electron at the gamma point of the
%Brillouin zone
%Putting it in units of cm^-3
Nc = 1/4.* (2.*meStar*me0.*kbT/(pi*hbar^2)).^(3/2)./(100)^3;

