function[Nc, Nv, ni, G, phi] =...
    calculateKeyParameters(temp, meStar, mhStar, photonEnergyum,...
nRefractive, kExtinction, eg300, eg,aGap,ncap,d, einf) 





%this has some similarity to Asbeck's calculation of photon recycling
%found here: P. Asbeck, J. Appl. Phy. vol 48, p 820 (1977). Self-absorption
%on the radiative lifetime in GaAs-GaAlAs double heterostructures
%doi: https://doi.org/10.1063/1.323633
photonRecycling = get_photon_recycling_factor(aGap, ncap, d);
phi = 1./(1-photonRecycling);

% Calculate the effective density of states in the valence and conduction
% band
Nc = densityInConduction(temp, meStar);
Nv = densityInValence(temp, mhStar);




%Calculate the intrinsic carrier concentration
ni = calculateIntrinsic(Nc, Nv,eg, temp);

%Calculate the radiative generation rate through the integration method
%%G = calcRadGenRate(nRefractive, kExtinction, temp,photonEnergyum, eg300);


% Calculate the radiative generation rate through an analytic form of the
% absorption coefficient using parabolic bands
G = calcAnalyticG(temp, meStar, mhStar, ni, eg,einf);
