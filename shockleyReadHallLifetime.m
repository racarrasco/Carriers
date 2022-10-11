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
% ni =          intrinsic carrier concentration (cm-3)
% Ev =          valence band edge (maximum) (eV)
% Ec =          conduction band edge (minimum) (eV)
% Et =          trap energy with respect to the conduction band edge (eV)
% sigmaN =      cross section multiplied by the trap density (m-1)
% dopingDensity doping density, can be either n-type or p-type units of
% (1e15 cm-3)
 

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
% Factor of 3 means that we are looking at the Root mean square of the 3D
% velocity the mean of the velocity has the 8/pi factor
% this will simply affect the defect concentration product, (for the most
% part) 
vp = sqrt(3*kb*x/(mhStar*meFree));
vn = sqrt(3*kb*x/(meStar*meFree));

%Lifetime in high majority carrier environment (low T as well) because
% minority carriers are low 
taup0 = 1./(sigmaN.*vp);
taun0 = 1./(sigmaN.*vn);




tauSRHRow = (taup0.*(n0 + n1) + taun0.*(p0 + p1) )./(p0 +n0);

%Make the lifetime in the end as a column vector
if(size(tauSRHRow,1)>1)
    tauSRH = tauSRHRow;
else
    tauSRH = tauSRHRow';
end