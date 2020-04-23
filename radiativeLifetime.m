function[tauRad] = radiativeLifetime(temp,a, b, c, d, e, f, g, hh,k)
%%%  Calculate the radiative lifetime  following the convention in a paper
%%%  by Olson: this is in the low excitation limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Input Parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% x         - The temperature of the sample
type = a; %p-doping or n-doping?
mestar = b;%electron effective mass in units of electron rest mass
mhstar = c;%hole effective mass in units of electron rest mass
eg = d; %band gap as a function of temperature
eg300 = e;%band gap at 300 K
refractive = f; %refractive index (n) of the material
extinction = g; %extinction coefficient (k) of the material
photonEnergy = hh; %photon energy in micron 
dopingDensity = k; %doping density in units of (cm)^-3

% k   =  The doping concentration in units of inverse cubic centimeters
if(size(temp,1)>1)
    x = temp';
else
    x = temp;
end
%effective electron and hole density of states as a function of temperature
Nc = densityInConduction(x, mestar); 
Nv = densityInValence(x, mhstar);

%Intrinsic electron density of states
ni = calculateIntrinsic(Nc,Nv,eg,x);
n0 = 0;
p0 = 0;
switch(type)
    case('p')
       %We have p-type doping
         %Uses the law of mass action to calculate the density of electrons
         % and holes in the conduction and valence band, respectively
         % (See: Steven H. Simon, The Oxford Solid State Basics
         
         p0 = 0.5.*(sqrt(dopingDensity^2 + 4.*ni.^2) + dopingDensity);
         n0 = ni.^2./p0;     
    case('n')
        %we have an n-type doping sample
        %uses the law of mass action to calculate the density of electrons
        %and holes in the conduction and valence band, respectively
        n0 = 0.5.*(sqrt(dopingDensity^2 + 4.*ni.^2) + dopingDensity);
        p0 = ni.^2./n0;
    otherwise
        %We have an intrinsic sample
        n0 = ni;
        p0 = ni;
        
        
end
    
G = calcRadGenRate(refractive, extinction, x, photonEnergy, eg300);
tauRadRow = ni.^2./(G.*(n0+p0));
if(size(tauRadRow,1)>1)
    tauRad = tauRadRow;
    
else
   tauRad = tauRadRow'; 
    
end

end