function[Nc, Nv, ni, G, phi] =...
    calculateKeyParameters(temp, meStar, mhStar, photonEnergyum,...
nRefractive, kExtinction, eg300, eg,aGap,ncap,d) 


photonRecycling = get_photon_recycling_factor(aGap, ncap, d);
phi = 1./(1-photonRecycling);
Nc = densityInConduction(temp, meStar);
Nv = densityInValence(temp, mhStar);
G = calcRadGenRate(nRefractive, kExtinction, temp,photonEnergyum, eg300);
ni = calculateIntrinsic(Nc, Np,eg, temp);
