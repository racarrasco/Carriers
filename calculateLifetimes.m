function[totalLifetime, tauRad, tauSRH, tauAug] = calculateLifetimes(...
    Tprobe, type, meStar, mhStar, egTempDep, eg300, valenceEdge, f1f2, conductionEdge,einf, nRefractive, kExtinction,...
    photonEnergyum,thickness, capRefractive, a, b, c)
% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, so a rename is performed in order to understand the code






dopingDensity = a;
defectLevel = b;
defectDensity = c;


%this has some similarity to Asbeck's calculation of photon recycling
%found here: P. Asbeck, J. Appl. Phy. vol 48, p 820 (1977). Self-absorption
%on the radiative lifetime in GaAs-GaAlAs double heterostructures
%doi: https://doi.org/10.1063/1.323633
photonRecycling = get_photon_recycling_factor(1200, capRefractive, thickness);

phi = 1/(1-photonRecycling);

%put it in units of microseconds
tauRad = phi.*1e6 .*radiativeLifetime(Tprobe, type, meStar, mhStar, egTempDep,...
    eg300, nRefractive, kExtinction, photonEnergyum, dopingDensity);

tauSRH = 1e6.*shockleyReadHallLifetime(Tprobe, type, meStar, mhStar,valenceEdge,...
                conductionEdge, defectLevel, defectDensity,dopingDensity);

tauAug = 1e6.*augerLifetime(Tprobe, type, meStar, mhStar, f1f2,egTempDep, einf,dopingDensity);

totalLifetime = 1./(1./tauRad + 1./tauSRH + 1./tauAug);
