function[radiativeGenCoefficient] = calcRadGenRate(n, k, temp, photonEnergyum, eg300)

% Calculate the radiative coefficient at a number of different temperatures
% The integral can be found in Blakemore's Semiconductor physics(1962)
% formula 511.3 

% Compare to eq 14 in W. V. Roosbroeck and W. Shockley, Phys. Rev. vol 94,
% 1558 (1954). Photon-Radiative Recombination of electron and Holes in Ge
% Doi: 10.1103/PhysRev.94.1558

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Input parameters %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%n                 =  Refractive index as a function of the wavelength
%k                 =  Extinction coefficient as a function of the
%wavelength
%temp              =  Input temperatures; this will be in the range of the
%measured lifetimes but with a finer mesh (unless a fit is being performed)

%photonEnergyum    =  Corresponding photon energy for the refractive
%index and extinction coefficient. The order of the photon energy must be
%in ascending order (decreasing energy)  and in units of microns

%eg300                = bandGap of the semiconductor at 300 K, this will be used to
%calculate the limits of integration to approximate the integral, Blakemore
%states that the integral will mostly matter at approximately the bandgap
%and then some 10kT above the bandGap

%sort in descending order for  ascending order in eV
[photonEnergyumSort, indices] = sort(photonEnergyum,'Descend');
nS = n(indices);
kS = k(indices);




%speed of light in cm/s
c = 2.99792e10;

% boltzmann constant in joules/Kelvin
boltzmann = 1.380649e-23;


%electron charge in Coulomb
e = 1.602176634e-19;  


%planck's constant in joules  second
h = 6.626070e-34;

%planck's constant in ev second
heV = h/e;

%Boltzmann constant in eV/Kelvin
boltzmanneV = boltzmann/e;

%put photon energy in units of cm
photonEnergycm = photonEnergyumSort.*1e-4;

%Photon energy in eV
photonEnergyeV = 1.24./photonEnergyumSort;

%absorption coefficient in inverse cm
alpha = 4*pi.*kS./photonEnergycm;

%thermal energy (kT)
kT = boltzmanneV.*temp;


%lower limit of integration is ~eg and higher limit is ~10kT higher than
%the band gap We are using the absorption coefficient at only 300K so the
%band gap is eg at ~300K
lowerEnergy = eg300 - 3*boltzmanneV*300;

%higherEnergy = eg300 + 15*boltzmanneV*300;

%indices of integration
intIndices = lowerEnergy<photonEnergyeV; %& ...
%    photonEnergyeV<higherEnergy;

%the integration values 
photonEnergyeVInteg = photonEnergyeV(intIndices);

%refractive index values
nInteg = nS(intIndices);

%absorption coefficient to be integrated
alphaInteg = alpha(intIndices);

%The denoninator for the integral, follows planck's harmonic oscillator for
%calculating the density of photons per unit volume
planck = exp(photonEnergyeVInteg./kT) - 1;


%Integrate to get the radiative generation rate
integrand = nInteg.^2.*alphaInteg.*photonEnergyeVInteg.^2./planck;
%figure()
%semilogy(photonEnergyeVInteg, integrand);

integResult = trapz(photonEnergyeVInteg, integrand);
radiativeGenCoefficient = 8*pi/heV^3/c^2.*integResult;

%Checking the units! (1/(ev s))^3  (s/cm)^2   (1/cm (ev)^2  (ev))
%                     1/h^3        (1/c)^2    (alpha  E^2    dE )
%                     factor       factor     integrand
%
% it checks out to 1/(s cm^3)   






