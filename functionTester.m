%% Lifetime Plotter
% Set up the parameters for testing

%Read in the band parameters
materialsDataBase = readtable('BandParameters.csv');
bandParams = materialsDataBase.InAs;

%band parameters and refractive index and absorption coefficient
%Set parameters for testing
[meStar, mhStar, eg, nRefractive,kExtinction, photonEnergyum, varshni, ...
     einf, f1f2] = extractParameters(bandParams, 'InAs');
 
%low doping
type = 'n';

%Create a higher temperature resolution
beginTemp = 0; %low limit temperature
endTemp = 900; %high limit temperature
number = endTemp - beginTemp +1; %steps of 1 Kelvin
Tprobe = linspace(beginTemp,endTemp, number);


eg300 = varshni(300); %band gap at 300 K
egTemp = varshni(Tprobe); %band gap as a function of temperature
conductionEdge = varshni(0); % conduction band edge is just band gap at 0 K

%Set the valence band edge as a 0 reference
valenceEdge = conductionEdge - egTemp;
dopingDensity = 3e16; %units of cm^(-3)

defectLevel = [5e-3, 25e-3, 40e-3]; %follow elizabeth's plots for SRH
defectDensity =70;  %use a defect density from Olson


%% Shockley read hall
figure()
for i = 1:length(defectLevel)
tauSRH = shockleyReadHallLifetime(Tprobe, type, meStar, mhStar, valenceEdge,...
          conductionEdge, defectLevel(i), defectDensity, dopingDensity);
name = string(int64(defectLevel(i)*1000)) + " meV";
semilogy(Tprobe,tauSRH,'DisplayName', name)
hold on
end
hold off
xlim([0,endTemp])
xlabel('Temperature (K)')
ylabel('Lifetime (s)')
lgd = legend;
lgd.Title.String = "E_T";
%ylim([1e-7, 1e-5])

title('SRH Lifetime: InAs')
%% Radiative recombination
dopingType = ["i", "n"];
name = ["intrinsic", "n-Type: " + string(dopingDensity) + " cm^{-3}"];
figure()
for i = 1:2
    
    tauRad = radiativeLifetime(Tprobe,dopingType(i),meStar, mhStar, eg,eg300,...
        nRefractive, kExtinction,photonEnergyum,dopingDensity);
    
    semilogy(Tprobe,tauRad,'Displayname',name(i))
    hold on
end
legend
xlabel('Temperature (K)')
ylabel('\tau_{rad} (s)')
ylim([1e-10, 1e-4])
xlim([beginTemp,endTemp])
title('InAs radiative lifetime')
%% Auger recombination 

figure()
for i = 1:2
    
    tauAug = augerLifetime(Tprobe, dopingType(i), meStar, mhStar, f1f2,egTemp,...
        einf,dopingDensity);
    
    semilogy(Tprobe,tauAug,'Displayname',name(i))
    hold on
end
legend
xlabel('Temperature (K)')
ylabel('\tau_{aug} (s)')
ylim([1e-10, 1])
xlim([beginTemp,endTemp])
title('InAs auger lifetime')

%% Putting it all together
tauall = 1./(1./tauRad + 1./tauSRH + 1./tauAug);
taus = [tauRad; tauSRH; tauAug; tauall];
names = ["Radiative"; "SRH"; "Auger";"total"];
figure()
for i =1:4
semilogy(Tprobe, taus(i,:), 'Displayname', names(i,:))
hold on 
end
legend
xlabel('Temperature (K)')
ylabel('Lifetime (s)')

