function[tauRad] = radiativeLifetime(type, ni, G, dopingDensity)
%%%  Calculate the radiative lifetime  following the convention in a paper
%%%  by Steenbergen: this is in the low excitation limit
%%%  TODO: will need to include higher excitation levels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Input Parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type     p-doping or n-doping?
% dopingDensity =  doping density in units of (cm)^-3
% G the radiative coefficient calculated in calcRadGenRate.m the G term is
% related to the B-Term in some literature, typically
% taurad = ni^2/G(n0 + po)   or
% taurad = 1/[B(n0 + p0)]

%Intrinsic electron density of states
n0 = 0;
p0 = 0;

% the input is in units of 1e15 1/cm^3 so lets put it in units of 1/cm^3
D = dopingDensity*1e15; 
switch(type)
    case('p')
       %We have p-type doping
         %Uses the law of mass action to calculate the density of electrons
         % and holes in the conduction and valence band, respectively
         % (See: Steven H. Simon, The Oxford Solid State Basics
         
         p0 = 0.5.*(sqrt(D^2 + 4.*ni.^2) + D);
         n0 = ni.^2./p0;     
    case('n')
        %we have an n-type doping sample
        %uses the law of mass action to calculate the density of electrons
        %and holes in the conduction and valence band, respectively
        n0 = 0.5.*(sqrt(D^2 + 4.*ni.^2) + D);
        p0 = ni.^2./n0;
    otherwise
        %We have an intrinsic sample
        n0 = ni;
        p0 = ni;
        
        
end

tauRadRow = ni.^2./(G.*(n0+p0));

%Make the lifetime a column vector
if(size(tauRadRow,1)>1)
    tauRad = tauRadRow;
    
else
   tauRad = tauRadRow'; 
    
end




end