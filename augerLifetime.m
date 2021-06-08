function[tauAuger] = augerLifetime(temp, a, b, c, d, f,g,h)
%Calculate the Auger process; it follows Chapter 6 in Blakemore, which in
%turn can be compared to reference: 
%B. V. Olson, E. A. Shaner, J. K. Kim, J. F. Klem, S. D. Hawkins, and M. E.
%Flatte, Appl. Phys. Lett. vol 103, 052106 (2013). Identification of
%dominant recombination mechanisms in narrow-bandgap InAs/InAsSb type-II
%superlattices and InAsSb alloys.
%doi: http://dx.doi.org/10.1063/1.4817400

%Original Auger paper: A. R. Beattie, and P. T. Landsberg, Proc. R. Soc.
%Lond. A. vol 249, pg 16-29 (1959).

%UPDATED ON 2/12/2021 TO INCLUDE HOLE COLLISIONS and it does not make any
%noticeable difference even with intrinsically p-type materials go ahead
%and try it for sample GN1886

%Make the temperature as a row because the band gaps are a row
if(size(temp,1) > 1)
    x = temp';
else
    x = temp; % the temperature is already a row
end


type = a; %The doping type
meStar = b; %electron effective mass
mhStar = c; %hole effective mass
F1F2 = d;  %Bloch overlap function
eg = f;    % temperature dependent band gap
einf = g;  %high-frequency dielectric constant]

%The original input is in units of 1e15 so lets put it to 1e15
dopingDensity = h*1e15; 

%electron charge in Coulomb
e = 1.602176634e-19;    


% boltzmann constant in joules/Kelvin
boltzmann = 1.380649e-23;

% boltzmann in ev/Kelvin
boltzmanneV = boltzmann/e;

%boltzmann energy
kT = boltzmanneV.*x;


NC = densityInConduction(x, meStar); %effective electron density of states
PC = densityInValence(x, mhStar); %effective hole density of states
ni = calculateIntrinsic(NC,PC,eg,x);% intrinsic carrier concentration



n0 = ni;
p0 = ni;


%Assume heavy holes and dominant mechanism are electron collisions (Auger-1) and pair
%creation from impact ionization 
mu = meStar/mhStar;
ratio = meStar;


switch(type)
    case('p')
       %We have p-type doping
         %Uses the law of mass action to calculate the density of electrons
         % and holes in the conduction and valence band, respectively
         p0 = 0.5.*(sqrt(dopingDensity^2 + 4.*ni.^2) + dopingDensity);
         n0 = ni.^2./p0;
    case('n')
        %We have n-type doping
        %uses the law of mass action to calculate the density of electrons
        %and holes in the conduction and valence band, respectively
        n0 = 0.5.*(sqrt(dopingDensity^2 + 4.*ni.^2) + dopingDensity);
        p0 = ni.^2./n0;
    otherwise
        %We have an intrinsic sample
        n0 = ni;
        p0 = ni;
        
        
end



%Determine the dominant Auger recombination type if the electrons are
%heavier, then the dominannt Auger is hole collisions (Auger 2)
%See chapter Chapter 6 in Blakemore, Semiconductor statistics
if(mu > 1)
    mu = mhStar/meStar;
    ratio  = mhStar;
end


tauAuger12 =  3.8e-18 * einf^2 * (1 + mu)^(1/2) * (1 + 2*mu) / (ratio * abs(F1F2)^2)...
    .*(eg./kT).^(3/2) .* exp((1 + 2*mu) / (1 + mu) .* eg./kT );

beta = mu^(1/2) * (1 + 2 * mu) / (2 + mu)*exp(- (1 - mu)/(1 + mu).*eg./kT);
tauAugerRow = 2*ni.^2 ./ ((n0 + p0).*(n0 + beta .* p0)) .*tauAuger12;
%%%%Old implementation
%%%%tauAugerRowee = 2*ni.^2 ./ ((n0 + p0).*(n0)) .*tauAuger12;


if(size(tauAugerRow,1)>1)
    %It's actually a column vector
    tauAuger = tauAugerRow;
else 
    tauAuger = tauAugerRow';
end

end


    
    
    
 
