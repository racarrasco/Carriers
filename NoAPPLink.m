%% Rather than using the app, just use this!!!
% This is for combining three different sample sets in order to minimize the
% amount of fit parameters

clear
%%% Change these parameters in order to change what data to fit.
samplesLoad = {'GN1881_FittedResults.mat';'GN1882_FittedResults.mat';...
    'GN1886_FittedResults.mat'};
plLoad = {'GN1881EMG_FittedResultsPL.mat';'GN1882EMG_FittedResultsPL.mat';...
    'GN1886EMG_FittedResultsPL.mat'};

sampleNames = cell(size(samplesLoad));
intentionalDoping = [false; true; true];
dopingTypes = ['n', 'n', 'p'];
dopingDensities = [1.0; 1.5; 2.5];
linkDopingLevel = [false; true; false];
f1f2 = 0.15;
meStar = 0.026; % meStar
mhStar = 0.3333333; % mhStar

conductionEdge = 0; % conduction Edge
einf = 12.2; % dielectric constant


defectLevelGuesses = [-0.067; -0.067; -0.067];
%defectLevel2Guesses = [NaN; -0.070;-0.070];
defectDensityGuesses = [4; 12; 4];
%defectDensity2Guesses = [NaN; 12; 4];


%doping in the samples
doping1 = dopingDensities(1);
doping2 = dopingDensities(2);
doping3 = dopingDensities(3);

%Defect levels in the samples
defectlevel1 = defectLevelGuesses(1);
defectlevel2 = defectLevelGuesses(2);
defectlevel3 = defectLevelGuesses(3);

%defect density in the samples
defectdens1 = defectDensityGuesses(1);
defectdens2 = defectDensityGuesses(2);
defectdens3 = defectDensityGuesses(3);




% Preallocate variables for loading the samples
% Array of excitations
excitationsArray = cell(size(samplesLoad));

% Array of lifetimes for the different samples
t0_vTempArray    = cell(size(samplesLoad));

% Array of temperatures for the different samples
temperaturesArray    = cell(size(samplesLoad));

% Array of thicknesses for the samples
thicknesses = [1000, 4000, 3500];

% Array of photon recycling factors
phis = zeros(size(samplesLoad));

% Array of data lengths
lengths = zeros(size(samplesLoad));

% Array of band gap lengths
varshniEinsteins = cell(size(samplesLoad));

% Load and store the data
for i = 1:size(samplesLoad,1)
    %Load lifetime data
    load(samplesLoad{i});
    %store excitations
    excitationsArray{i} = Excitations;
    %store lifetime results from lowest excitation
    t0_vTempArray{i} = t0_vTemp(1,:,end);
    
    sampleNames{i} = SampleName;
    
    %store the temperature
    temperaturesArray{i} = Temperatures; 
    
    %lengths of temperatures
    lengths(i) = size(Temperatures,2);
    
    %photon recycling factors
    phis(i) = 1./(1-get_photon_recycling_factor(2550, 3.7,thicknesses(i)));
    
    %load PL results
    load(plLoad{i});
    varshniEinsteins{i} = @(T) boseEinstein(T, boseParameters(1,1),...
        boseParameters(2,1), boseParameters(3,1));
end
allTemps = zeros(sum(lengths),1);
allLTs = zeros(sum(lengths),1);
bandGaps = zeros(sum(lengths),1);





%concatenate the temperature data and lifetime data
for i = 1:size(lengths,1)
    if(i == 1)
       allTemps(1:lengths(i)) = temperaturesArray{i};
       allLTs(1:lengths(i)) = t0_vTempArray{i};
       bandGaps(1:lengths(i)) = varshniEinsteins{i}(temperaturesArray{i})./1000;
    else
        lastIndex = sum(lengths(1:i - 1));
        allTemps(lastIndex + 1:lastIndex +lengths(i)) = temperaturesArray{i};
        allLTs(lastIndex +1:lastIndex + lengths(i)) = t0_vTempArray{i};
        bandGaps(lastIndex + 1:lastIndex + lengths(i)) = ...
        varshniEinsteins{i}(temperaturesArray{i})./1000;
    end
end

% valence band edge
valenceEdge = conductionEdge - bandGaps;

%{
totalLifetime = calculateTotalLifetimeLink(...
    tProbeX, type, bandGaps, valenceEdge, conductionEdge, einf, phi(1), phi(2),...
    meStar, mhStar, f1f2, dopingDensity1, dopingDensity2, defectLevelGuess,...
    defectDensity1, defectDensity2, lengthOne);
%}


%{  
%HERE
totalLifetimeLink = calculateTotalLifetimeLink(allTemps, lengths,...
    dopingTypes, bandGaps', valenceEdge', 0, einf, phis,...
    meStar, mhStar, f1f2, doping1, doping2, doping3,...
    defectlevel1,...



%{    
defectlevel2,
 %}
    defectlevel3,...
    defectdens1, defectdens2, defectdens3);
AND HERE
%}

%% Fit with coupled doping

[fobjects, gofs, outputs] = fitLifetimeLink(allTemps, allLTs,lengths, dopingTypes, bandGaps, valenceEdge, conductionEdge,...
    einf, phis, meStar, mhStar, f1f2, doping1, doping2, doping3, defectlevel1,defectlevel3,...
    defectdens1, defectdens2, defectdens3);
%%
successfulIndices = [];
successLength = 0;
for i =1:size(outputs,1)
    messagei = outputs(i).message;
    if contains(messagei, 'Success')
        successLength = successLength + 1;
     successfulIndices(successLength) = i;
    
    end

end




%% Fit with decoupled doping

[fobjects2, gofs2, outputs2] = fitLifetimeLink(allTemps, allLTs,lengths, dopingTypes, bandGaps, valenceEdge, conductionEdge,...
    einf, phis, meStar, mhStar, f1f2, doping1, doping2, doping3, defectlevel1,defectlevel3,...
    defectdens1, defectdens2, defectdens3, defectlevel2);
%% Fit with fixed effective masses






%% Plot all of the data in subplots
xin = allTemps;
yin = allLTs;



iterations = size(fobjects,1);
row = 1;
column = 1;
numSubPlots = 25;
for i  = 1:iterations
    %make a 5X5 subplot where each plot is the best fit
    remainder = rem(i, numSubPlots);
    
    %The remainder is the the subplot index
    if (remainder == 1)
        figure
        %calculate whats left
        leftToGo = iterations - i;
        
        if(leftToGo ~= 0)
            column = 5;
            row = 5;
        else
            column = 1;
            row  = 1;
        end
    
    %no remainder, then index the last subplot
    elseif(remainder==0)
        
        remainder = numSubPlots;
        
    end
    
    %the goodness of fit value
    gofi = gofs(i);
    %subplot
    subplot(row,column,remainder)
    
    fobjecti = fobjects{i};
    
    yi = fobjecti(xin);

% show the data
 plot(xin(1:lengths(1)), yin(1:lengths(1)),'o','DisplayName','undoped')
    hold on
    plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
        yin(lengths(1) + 1 :sum(lengths(1:2))),'+','DisplayName','n-type')
    lastIndex = sum(lengths(1:2));
    
    plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
        yin(lastIndex + 1: lastIndex+lengths(3)),'>','DisplayName','p-type');
    
    
    %Show the fits
    plot(xin(1:lengths(1)), yi(1:lengths(1)),'k-','DisplayName','fit')
    
    plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
        yi(lengths(1) + 1 :sum(lengths(1:2))),'k--','DisplayName','fit')
    
    plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
        yi(lastIndex + 1: lastIndex+lengths(3)),'k.','DisplayName','fit');
    titleString = "i = " +string(i) +"/" +string(iterations) + "; rmse: " + string(gofi.rmse);
    title(titleString);
%{
    xlabel("Temperature (K)");
    ylabel("Minority carrier lifetime (us)");
  %} 
    
    hold off

end


%% PLot all lifetime mechanisms to all!


numPlots = size(fobjects,1);
numSubPlots = 25;
row = 1;
column = 1;
[~, indBest] = min([gofs.rmse]);
for indID = 1:numPlots
    
    
     %make a 5X5 subplot where each plot is the best fit
    remainder = rem(indID, numSubPlots);
    
    %The remainder is the the subplot index
    if (remainder == 1)
        figure
        %calculate whats left
        leftToGo = numPlots - indID;
        
        if(leftToGo ~= 0)
            column = 5;
            row = 5;
        else
            column = 1;
            row  = 1;
        end
    
    %no remainder, then index the last subplot
    elseif(remainder==0)
        
        remainder = numSubPlots;
        
    end
    
    %the goodness of fit value
    gofi = gofs(indID).rmse;
    %subplot
    subplot(row,column,remainder)
    
    

% show the data
Tprobe = linspace(50,330,281);
fobjecti = fobjects{indID};
reCalcegTempDep = varshniEinstein(Tprobe)./1000;
reCalcNc = densityInConduction(fobjecti.b,Tprobe);
reCalcNv = densityInValence(fobjecti.c, Tprobe);
reCalcni = calculateIntrinsic(reCalcNc,reCalcNv, reCalcegTempDep, Tprobe);
reCalcG  = calcAnalyticG(Tprobe, fobjecti.b, fobjecti.c, reCalcni, reCalcegTempDep,...
    fobjecti.g);

if(exist('tauSRH2'))
    [calctot, calcRadm, calcSRH, calcSRH2, calcAug] = ...
    calculateLifetimes2(Tprobe, fobjecti.a,fobjecti.b, fobjecti.c, reCalcegTempDep,...
    0-reCalcegTempDep,0, fobjecti.g, reCalcNc,...
    reCalcNv,reCalcni,reCalcG , phi,fobjecti.o,fobjecti.p, fobjecti.q,fobjecti.r,...
    fobjecti.s, fobjecti.t);
else
[calctot, calcRadm, calcSRH, calcAug] = ...
    calculateLifetimes(Tprobe, fobjecti.a,fobjecti.b, fobjecti.c, reCalcegTempDep,...
    0-reCalcegTempDep,0, fobjecti.g, reCalcNc,...
    reCalcNv,reCalcni,reCalcG , phi,fobjecti.o,fobjecti.p, fobjecti.q,fobjecti.r);
end

semilogy(Temperatures, t0_vTemp(1,:,end),'o','Displayname', 'Data')
hold on
plot(Tprobe, calctot,'-', 'Displayname', 'Total')
plot(Tprobe, calcSRH,'--','Displayname','SRH')
plot(Tprobe, calcAug, '-.','Displayname','Auger')
plot(Tprobe, calcRadm, '--','Displayname','Radiative')
hold off
titlename = strcat('n0 = ',num2str(fobjecti.p),';',...
    'Et = ', num2str(-1000*fobjecti.q),...
    ' rmse = ',num2str(gofi));
if(indID == indBest)
   titlename = strcat(titlename,',','Best'); 
    
end
title(titlename)
xlim([70,300])
ylim([0.05,200])

  
%{
    xlabel("Temperature (K)");
    ylabel("Minority carrier lifetime (us)");
  %} 
    
  
    
    
end






%% Choose and blow up a single subplot

% This is a copy paste of what was in fitLifetime link
xin = allTemps;
yin = allLTs;
iterations = size(theObjects,1);

gofi = theGofs(index);
    
fobjecti = theObjects{index};

yi = fobjecti(xin);
  
figure

%%% show the data %%%
% plot the first section of the data which is the undoped data set
plot(xin(1:lengths(1)), yin(1:lengths(1)),'o','DisplayName','undoped')
hold on
% plot the n-type data
plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
     yin(lengths(1) + 1 :sum(lengths(1:2))),'+','DisplayName','n-type')

%The last data set
lastIndex = sum(lengths(1:2));

% plot the p-type data
plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
    yin(lastIndex + 1: lastIndex+lengths(3)),'>','DisplayName','p-type');pea
    
%Show the fits for the undoped data
plot(xin(1:lengths(1)), yi(1:lengths(1)),'k-','DisplayName','fit')

%Show the fit for the n-type data
plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
     yi(lengths(1) + 1 :sum(lengths(1:2))),'k--','DisplayName','fit')
 
%Show the fit for the p-type data    
plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
     yi(lastIndex + 1: lastIndex+lengths(3)),'k.','DisplayName','fit');

titleString = "i = " +string(index) +"/" +string(iterations) + "; rmse: " + string(gofi.rmse);
title(titleString);