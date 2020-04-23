function[ mestar, mhstar, eg, nRefractive, kExtinction, photonEnergyum, varshni, einf, f1f2] = extractParameters(bandParams, material)

pathAndFile = "OpticalConstantsDataBase\" + string(material) +".txt";


 %Open the .txt file yourself and see how the data is formatted!!!
 %We only have optical constants for InAs and InSb at room temperature (Not good)
    %It would be better to create temperature dependent optical constants
switch(material)
    case('InAs')
        opticalConstantsTable = readtable(pathAndFile);
        opticalConstantsTable.Properties.VariableNames = {'um', 'n', 'k'};
    case('InSb')
        opticalConstantsTable = readtable(pathAndFile);
        opticalConstantsTable.Properties.VariableNames = {'eV','n', 'k'};
        opticalConstantsTable.um = 1.24./opticalConstantsTable.eV;
        
        %USE InAsSb optical constants, but use InAs effective mass and
        %other band parameters
    case('InAsSb')
        opticalConstantsTable = readtable(pathAndFile,'HeaderLines', 3);
        opticalConstantsTable.Properties.VariableNames = {'um', 'n','k'};
    case('GaAs')
        opticalConstantsTable = readtable(pathAndFile);
        opticalConstantsTable.Properties.VariableNames = {'eV','e1', 'e2'};
        opticalConstantsTable.um = 1.24./opticalConstantsTable.eV;
        e1 = opticalConstantsTable.e1;
        e2 = opticalConstantsTable.e2;
        opticalConstantsTable.n = 1/sqrt(2) .*sqrt(e1 + sqrt(e1.^2 +e2.^2));
        opticalConstantsTable.k = 1/sqrt(2).*sqrt(-e1 + sqrt(e1.^2 + e2.^2));
end

photonEnergyum = opticalConstantsTable.um;
kExtinction = opticalConstantsTable.k;
nRefractive = opticalConstantsTable.n;

%row1     =  Band gap in eV (0 K) perhaps
%row2     =  electron effective mass in units of free electron mass
%row3     =  gamma1 inverse effective mass parameter 1
%row4     =  gamma2 inverse effective mass parameter 2
%row5     =  gamma3 inverse effective mass parameter 3
%row6     =  Varshni equation parameter1
%row7     =  Varshni equation parameter2
%row8     =  high frequency dielectric constant
%row9     =  Bloch overlap integral (not as a function of temperature)

eg = bandParams(1);
mestar  = bandParams(2);
gamm1 = bandParams(3);
gamm2 = bandParams(4);

%The heavy hole effective mass is a function of the inverse effective mass
%parameters The formula for calculating the heavy hole effective mass is in
% equation 2.16 in Vurgaftman's paper J. Appl. Phys. vol 89, 5815 (2001)
%Band Parameters for III-V compound semiconductors and their alloys
mhstar = 1/(gamm1 - 2*gamm2);

varshni = @(T) eg - bandParams(6).*T.^2./(T + bandParams(7));
einf = bandParams(8);
f1f2 = bandParams(9);


     
        


