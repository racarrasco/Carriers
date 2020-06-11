function[Nc, Nv, ni, G] = calculateKeyParameters2(temp, meStar, mhStar,...
    eg, einf) 

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
