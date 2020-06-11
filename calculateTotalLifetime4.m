function[totalLifetime] = calculateTotalLifetime4(...
    x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)
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
% b = egTempDep the temperature dependent band gap
% c = Ev the valence edge energy (taken as the negative band gap)
% d = Ec the conduction edge energy (taken as the 0 reference) 
% e  = einf the static dielectric constant below the absorption edge
% f  = phi the photon recycling factor


% FITTING PARAMETERS
% g = mestar 
% h = mhstar
% k = Auger overlap parameter  OPTIONAL
% l = doping density in units of 1e-15cm^(-3)
% m = defect level in reference to the conduction band
% n = defect density in units of inverse meters

% o = defect level 2 in reference to the conduction band
% p =  defect density 2  in units of inverse meters

[Nc, Nv, ni, G] = calculateKeyParameters2(x, g, h,b, e);




%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, so a rename is performed in order understand the code
dopingDensity = l;
defectLevel = m;
defectDensity = n; 




%put it in units of microseconds

%tauRad = phi(f).*1e6 .*radiativeLifetime(type (a), ni, G, dopingDensity);
tauRad = f.*1e6 .*radiativeLifetime(a, ni, G, dopingDensity);


% tauSRH =  
% 1e6.*shockleyReadHallLifetime(Tprobe, type(a), meStar(g), mhStar(h),Nc (h),
%                Nv, ni, valenceEdge(c), conductionEdge(d), defectLevel, 
%                defectDensity, dopingDensity);
tauSRH = 1e6.*...
       shockleyReadHallLifetime(x, a, g, h, Nc, Nv, ni, c, d, defectLevel,...
       defectDensity, dopingDensity);

tauSRH2 = 1e6.*...
       shockleyReadHallLifetime(x, a, g, h, Nc, Nv, ni, c, d, o,...
       p, dopingDensity);


% tauAug = 1e6.*augerLifetime(Tprobe(x), type(a), meStar(g), mhStar(h), f1f2(o),
% egTempDep(b), einf(e),dopingDensity);
tauAug = 1e6.*augerLifetime(x, a, g, h, k, b, e, dopingDensity);



totalLifetime = 1./(1./tauRad + 1./tauSRH + 1./tauAug + 1./tauSRH2);
