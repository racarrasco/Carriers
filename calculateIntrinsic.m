function[intrinsic] = calculateIntrinsic(Nc, Np,eg, T)
%Calculate the intrinsic carrier concentration as a function of temperature

% boltzmann constant in joules/Kelvin
boltzmann = 1.38064852e-23;

%electron charge in Coulomb
e = 1.602176634e-19;  

%Boltzmann constant in eV/Kelvin
boltzmanneV = boltzmann/e;

kT = boltzmanneV.*T;

%In units of cm^(-3)  because Nc and Np are already in units of cm^(-3)
intrinsic = sqrt(Nc.*Np.*exp(-eg./kT));