function[tauSRH] = ...
    shockleyReadHallLifetime(temp, type,meStar,mhStar, Nc, Nv, ni...
    ,Ev,Ec,Et,sigmaN,dopingDensity)

% Calculate the Shockley-Read hall lifetime as a function of temperature

% See equation 5.5 from
% W. Shockley and W. T. Read Jr, Phys. Rev. Vol 87, p. 835 (1952).
% Statistics of the recombination of Holes and Electrons
% Doi: https://doi.org/10.1103/PhysRev.87.835
% The equation assumes a low excess carrier density


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Input parameters %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type =        p-type doping or n-type doping
% meStar =      electron effective mass
% mhStar =      hole effective mass
% Nc =          effective density of States in the conduction band
% Nv =          effective density of States in the valence band
% ni =          intrinsic carrier concentration
% Ev =          valence band edge (maximum)
% Ec =          conduction band edge (minimum)
% Et =          trap energy
% sigmaN =      cross section multiplied by the trap density
% dopingDensity doping density, can be either n-type or p-type
 

ne = 1e15;
%make the temperature a row since the band gaps are a row-wise vector
if(size(temp,1) > 1)
    x = temp';
else
    x = temp;
end


%Boltzmann constant in Joules/Kelvin
kb = 1.380649e-23;

%electron charge in Coulomb
e = 1.602176634e-19;

%Boltzmann constant in eV/Kelvin
kbeV = kb/e;


%BoltzmanEnergy in eV (row -wise vector)
kbTeV = kbeV.* x;

%Free electron mass in kg
meFree = 9.10938356e-31;

n0 = 0;
p0 = 0;

%the input is in units of 1e15 1/cm^3, so we need to put it in 1/cm^3
D = dopingDensity*1e15;
switch(type)
    case('p')
       %We have p-type doping
         %Uses the law of mass action to calculate the density of electrons
         %and holes in the conduction and valence band, respectively
         p0 = 0.5.*(sqrt(D^2 + 4.*ni.^2) + D);
         n0 = ni.^2./p0;
    case('n')
        %we have n-type doping
        %uses the law of mass action to calculate the density of electrons
        %and holes in the conduction and valence band, respectively
        n0 = 0.5.*(sqrt(D^2 + 4.*ni.^2) + D);
        p0 = ni.^2./n0;
    otherwise
        %We have an intrinsic sample
        n0 = ni;
        p0 = ni;        
end


% Electron and hole density in the conduction and valence bands,
% respectively (row-wise vector)
n1 = Nc.*exp((Et - Ec)./kbTeV);
p1 = Nv.*exp((Ev - Et)./kbTeV);

% Thermal velocities of electron and holes
% formula Taken from
% E. H. Steenbergen, B. C. Connelly, G. D. Metcalfe, H. Shen, M. Wraback,
% D. Lubyshev, Y. Qiu, J. M. Fastenau, A. W. K. Liu, S. Elhamri, O. O.
% Cellek, and Y.-H. Zhang, Proc. SPIE 8512, 85120L (2012)
% DOI: 10.1117/12.930949tau
vp = sqrt(8*kb*x/(pi*(mhStar*meFree)));
vn = sqrt(8*kb*x/(pi*(meStar*meFree)));


taup0 = 1./(sigmaN.*vp);
taun0 = 1./(sigmaN.*vn);




tauSRHRow = (taup0.*(n0 + n1 + ne) + taun0.*(p0 + p1+ne) )./(p0 +n0 + ne);

%Make the lifetime in the end as a column vector
if(size(tauSRHRow,1)>1)
    tauSRH = tauSRHRow;
else
    tauSRH = tauSRHRow';
end