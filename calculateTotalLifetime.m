function[totalLifetime] = calculateTotalLifetime(...
    x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s)
% Calculate all of the lifetime components and add them together. The rates
% are added then the reciprocal of the total rate is the total lifetime
% This is for the fitting routine where the single output is the total
% lifetime output, SO THIS FILE WILL BE HARD TO READ, IN ORDER TO
% UNDERSTAND THE CODE BETTER (AND QUICKER) SEE THE FILE 
% calculateLifetimes.m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = Temperature that will be used to calculate the total lifetime
% a = doping type n-type or p-type
% b = meStar, the effective mass of conduction electron
% c = mhStar, the effective mass of the heavy hole valence electron
% d = egTempDep the temperature dependent band gap
% e = Ev the valence edge energy (taken as the negative band gap)
% f = Ec the conduction edge energy (taken as the 0 reference) 
%%%% g = f1f2 the overlap integral for auger calculation
%g  = einf the static dielectric constant below the absorption edge
%h  = Nc the conduction band effective density of states
%k  = Nv the valence band effective density of states
%l  = ni the intrinsic carrier concentration
%m  = G  the radiative generation rate 
%n  = phi the photon recycling factor


% FITTING PARAMETERS
% o = Auger overlap parameter  OPTIONAL
% p = doping density in units of 1e-15cm^(-3)
% q = defect level in reference to the conduction band (eV)
% r = defect density in units of inverse meters
% s = cross section activation energy (eV)



%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, so a rename is performed in order understand the code
dopingDensity = p;
defectLevel = q;
defectDensity = r; 
crossSection1  = s;




%put it in units of microseconds

%tauRad = phi(n).*1e6 .*radiativeLifetime(type (a), ni(l), G(m), dopingDensity);
tauRad = n.*1e6 .*radiativeLifetime(a, l, m, dopingDensity);


% tauSRH =  
% 1e6.*shockleyReadHallLifetime(Tprobe, type(a), meStar(b), mhStar(c),Nc (h),
%                Nv(k), ni(l), valenceEdge(e), conductionEdge(f), defectLevel, 
%                defectDensity, dopingDensity);
tauSRH = 1e6.*...
       shockleyReadHallLifetime(x, a, b, c, h, k, l, e, f, defectLevel,...
       defectDensity, dopingDensity, crossSection1);


% tauAug = 1e6.*augerLifetime(Tprobe(x), type(a), meStar(b), mhStar(c), f1f2(o),
% egTempDep(d), einf(g),dopingDensity);
tauAug = 1e6.*augerLifetime(x, a, b, c, o, d, g, dopingDensity);



totalLifetime = 1./(1./tauRad + 1./tauSRH + 1./tauAug);
