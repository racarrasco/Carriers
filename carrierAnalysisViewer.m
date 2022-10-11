%% Lifetime Plotter
% Set up the parameters for testing

%Read in the band parameters
tableDataBase = 'BandParameters.csv';
materialsDataBase = readtable(tableDataBase);
bandParams = materialsDataBase.InAs;

%band parameters and refractive index and absorption coefficient
%Set parameters for testing
%[meStar, mhStar, eg, nRefractive,kExtinction, photonEnergyum, varshniEinstein, ...
 %    einf, f1f2] = extractParameters(bandParams, 'InAs');
meStar = 0.026;
mhStar = 0.3333;

%low doping
type = 'p';

%Create a higher temperature resolution
beginTemp = 50; %low limit temperature
endTemp = 330; %high limit temperaturef
number = endTemp - beginTemp +1; %steps of 1 Kelvin
Tprobe2 = linspace(beginTemp,endTemp, number);


eg300 = varshniEinstein(300)/1000; %band gap at 300 K
egTemp = varshniEinstein(Tprobe2)/1000; %band gap as a function of temperature
conductionEdge = 0; % conduction band edge is 0 energy reference

%Set the valence band edge as negative band gap
valenceEdge = conductionEdge - egTemp;
dopingDensity = [bestInputs.Values(2), 3]; %units of x10^(15)cm^(-3)

defectLevel = [bestInputs.Values(3), bestInputs.Values(3)]; 
defectDensity = [bestInputs.Values(4), 10];



%% Shockley read hall
figure()
tauSRHtest =  zeros([length(Tprobe2), length(defectLevel)]);
Nc = densityInConduction(Tprobe2, meStar);
Nv = densityInValence(Tprobe2, mhStar);
ni = calculateIntrinsic(Nc, Nv,egTemp, Tprobe2);
for i = 1:length(defectLevel)
tauSRHtest(:,i) = 1e6.*shockleyReadHallLifetime(Tprobe2, type, meStar, mhStar,Nc,Nv,ni,...
    valenceEdge, conductionEdge, defectLevel(i), defectDensity(i), dopingDensity(i));
name = sprintf("%0.2f meV", defectLevel(i)*1000);
plot(Tprobe2,tauSRHtest(:,i),'-','DisplayName', name)
hold on
end
plot(Tprobe, totalLifetime,'k-', 'DisplayName', 'Lifetime fit')
plot(Temperatures, t0_vTemp(1,:,end),'o','DisplayName','Data')
hold off
xlim([70,endTemp])
xlabel('Temperature (K)')
ylabel('Lifetime (s)')
lgd = legend;
%lgd.Title.String = "E_T";
%ylim([1e-7, 1e-5])

title('SRH Lifetime: InAs')
%% Radiative recombination
dopingType = ["i", "n"];
name = ["intrinsic", "n-Type: " + string(dopingDensity) + " cm^{-3}"];
figure()
for i = 1:2
    
    tauRad = radiativeLifetime(Tprobe,dopingType(i),meStar, mhStar, eg,eg300,...
        nRefractive, kExtinction,photonEnergyum,dopingDensity);
    
    plot(Tprobe,tauRad,'Displayname',name(i))
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
tauall = 1./(1./tauRad + 1./tauSRHtest + 1./tauAug);
taus = [tauRad; tauSRHtest; tauAug; tauall];
names = ["Radiative"; "SRH"; "Auger";"total"];
figure()
for i =1:4
semilogy(Tprobe, taus(i,:), 'Displayname', names(i,:))
hold on 
end
legend
xlabel('Temperature (K)')
ylabel('Lifetime (s)')

