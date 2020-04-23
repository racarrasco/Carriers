function[totalLifetime] = calculateTotalLifetime(...
    x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s)
% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% This is for the fitting routine where the single output is the total
% lifetime output, This is 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, so a rename is performed in order understand the code
dopingDensity = q;
defectLevel = r;
defectDensity = s; 


%this has some similarity to Asbeck's calculation of photon recycling
%found here: P. Asbeck, J. Appl. Phy. vol 48, p 820 (1977). Self-absorption
%on the radiative lifetime in GaAs-GaAlAs double heterostructures
%doi: https://doi.org/10.1063/1.323633
photonRecycling = get_photon_recycling_factor(1200, p, o);

phi = 1/(1-photonRecycling);

%put it in units of microseconds
tauRad = phi.*1e6 .*radiativeLifetime(x, a, b, c, d,...
    e, l, m, n, dopingDensity);

tauSRH = 1e6.*shockleyReadHallLifetime(x, a, b, c,f,...
                h, defectLevel, defectDensity,dopingDensity);

tauAug = 1e6.*augerLifetime(x, a, b, c, g,d, k,dopingDensity);

totalLifetime = 1./(1./tauRad + 1./tauSRH + 1./tauAug);
