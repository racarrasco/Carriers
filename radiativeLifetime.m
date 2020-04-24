function[tauRad] = radiativeLifetime(temp,type, mestar, mhstar,ni, G, eg300, k)
%%%  Calculate the radiative lifetime  following the convention in a paper
%%%  by Olson: this is in the low excitation limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Input Parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% x         - The temperature of the sample
%type     p-doping or n-doping?
mestar = b;%electron effective mass in units of electron rest mass
mhstar = c;%hole effective mass in units of electron rest mass
eg = d; %band gap as a function of temperature

dopingDensity = k; %doping density in units of (cm)^-3

% k   =  The doping concentration in units of inverse cubic centimeters
if(size(temp,1)>1)
    x = temp';
else
    x = temp;
end

%Intrinsic electron density of states
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
    
tauRadRow = ni.^2./(G.*(n0+p0));



if(size(tauRadRow,1)>1)
    tauRad = tauRadRow;
    
else
   tauRad = tauRadRow'; 
    
end

end