function[ mestar, mhstar, eg, nRefractive, kExtinction, photonEnergyum, varshni, einf, f1f2, ag]... 
    = extractParameters(bandParams, material)
%Extract the band parameters necessary to calculate the lifetime components

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Input Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bandParams = an empt



%Point to the file with 
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
%row3     =  gamma1 inverse effective mass parameter 1 Vurgaftman review
%row4     =  gamma2 inverse effective mass parameter 2 Vurgaftman review
%row5     =  gamma3 inverse effective mass parameter 3 Vurgaftman review

%row6     =  Varshni equation parameter1 Vurgaftman review
%row7     =  Varshni equation parameter2 Vurgaftman review
%row8     =  high frequency dielectric constant Zollner JVSTB 2019
%row9     =  Bloch overlap integral (not as a function of temperature) a


%row10    =  Absorption coefficient of the material at the band gap energy
%            Comes from S. T. Schaefer, S. Gao, P. T. Webster, R. R.
%            Kosireddy, and S. R. Johnson, J. Appl. Phys. vol 27 p 165705
%            (2020). Absorption edge characteristics of GaAs, GaSb, InAs,
%            InSb

eg = bandParams(1);
mestar  = bandParams(2);
gamm1 = bandParams(3);
gamm2 = bandParams(4);

%The heavy hole effective mass is a function of the inverse effective mass
%parameters The formula for calculating the heavy hole effective mass is in
% equation 2.16 in Vurgaftman's paper J. Appl. Phys. vol 89, 5815 (2001)
%Band Parameters for III-V compound semiconductors and their alloys
multiplier = 1;

mhstar = multiplier.*1/(gamm1 - 2*gamm2);

varshni = @(T) eg - bandParams(6).*T.^2./(T + bandParams(7));
einf = bandParams(8);
f1f2 = bandParams(9);
ag = bandParams(10);


     
        


