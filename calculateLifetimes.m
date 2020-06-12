function[totalLifetime, tauRad, tauSRH, tauAug] =...
    calculateLifetimes(Tprobe, type, meStar, mhStar, egTempDep, ...
    valenceEdge, conductionEdge,einf,Nc, Nv, ni, G, phi, f1f2, a, b, c,d)
% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These will be the fit parameters, and MATLAB arranges the fit parameters
% alphabetically, so a rename is performed in order to understand the code
dopingDensity = a;
defectLevel = b;
defectDensity = c;
activationEnergy = d;

% this has some similarity to Asbeck's calculation of photon recycling
% found here: P. Asbeck, J. Appl. Phy. vol 48, p 820 (1977). Self-absorption
% on the radiative lifetime in GaAs-GaAlAs double heterostructures
% doi: https://doi.org/10.1063/1.323633

% put it in units of microseconds
% radiative lifetime2
tauRad = phi.*1e6 .*radiativeLifetime(type, ni, G, dopingDensity);

% Shockley-Read-Hall lifetime
tauSRH = ...
1e6.*shockleyReadHallLifetime(Tprobe, type, meStar, mhStar,Nc, Nv, ni,...
valenceEdge, conductionEdge, defectLevel, defectDensity,dopingDensity,...
activationEnergy);

% Auger recombination lifetime
tauAug = 1e6.*augerLifetime(Tprobe, type, meStar, mhStar, f1f2,egTempDep,...
    einf,dopingDensity);


% Total lifetime as a function of temperature
totalLifetime = 1./(1./tauRad + 1./tauSRH + 1./tauAug);
