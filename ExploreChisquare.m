%% Look at deviations of the chi-square
filenameB = 'GN1881-combined11102020_LowestExcitation-03-16-2021-correctweights_CarrierAnalysis';
filenameLtsB = 'GN1881-combined11102020_FittedResultsLT';
load(filenameB,'bestInputs','indexBest','fobjects')
load(filenameLtsB,'t0_vTemp','t0Sigma_vTemp','Temperatures')
% Doping in units of 1e15
dopelow = 0.01;
dopehigh = 0.3;
% Defect level in units of eV below the conduction band  which is the reference
defectlow = -0.115;
defecthigh = -0.085;
% Resolution of the fits
resDope = 60;
resDefect = 30;
%{
[GN1881fobjects,GN1881Doping, GN1881Defect,GN1881DChiSquare] = lookAtChiSquare(filenameB,...
    filenameLtsB,...
    dopelow,dopehigh,defectlow,defecthigh,resDope, resDefect);
%}
%% Take 1 sigma and 2 sigma from the doping value (For GN1881)
% 1 sigma and then 2 sigma repeated
dopingInterests = [bestInputs.Values(2) - bestInputs.Errors(2),...
    bestInputs.Values(2) + 2*bestInputs.Errors(2),...
    bestInputs.Values(2) + 2*bestInputs.Errors(2)];
% doping densities from low doping to high doping range
dopingDensities = linspace(dopelow,dopehigh, resDope);

%The doping spacing
dopingSpacing = (dopehigh-dopelow)/resDope;

%set an n-by-length(dopingInterest) to designate the index where we have
%the 1 sigma and 2 sigma from best fit
dopedIndices = true([length(dopingDensities),length(dopingInterests)]);

%find the index that is closet to our one and two sigma
for i = 1:length(dopingInterests)
dopedIndices(:,i)= dopingInterests(i)-dopingSpacing*0.5 <=dopingDensities &...
    dopingDensities <= dopingInterests(i)+ dopingSpacing*0.5;
end

%Take one sigma and two sigma for the defect level and then 90 meV below
%the conduction band 
defectInterests =[bestInputs.Values(3) - bestInputs.Errors(3),...
    bestInputs.Values(3) + 2*bestInputs.Errors(3),...
    -0.090];
% Defect levels used in contour
defectLevels = linspace(defectlow,defecthigh,resDefect);

% Spacing
defectSpacing = (defecthigh - defectlow)/resDefect;

% Indices for the defect column 1 will be 1 sigma, then 2 sigma, then 90
% meV below conduction band
defectIndices = true([length(defectLevels),length(defectInterests)]);

% Find the indices closest to defectInterests
for i = 1:length(defectInterests)
defectIndices(:,i)= defectInterests(i)- defectSpacing*0.5 <=defectLevels &...
    defectLevels <= defectInterests(i) + defectSpacing*0.52;
end
%Instantiate fit objects with the n-sigma from doping and defect levels
fobjectInterest = cell(length(defectInterests),1);

for i =1: length(fobjectInterest)
    subplot(2,2,i + 1)
    fobjectInterest{i} =  GN1881fobjects{defectIndices(:,i), dopedIndices(:,i)};
    %figure('Name',strcat('Index Number:',num2str(i)))
    errorbar(Temperatures, t0_vTemp(1,:,end),t0Sigma_vTemp(1,:,end),'o','Displayname','Data')
    hold on
    plot(Temperatures,fobjects{indexBest}(Temperatures),'DisplayName','Best Fit (in paper)')
    plot(Temperatures, fobjectInterest{i}(Temperatures),...
        'DisplayName','newFit')
    hold off
    legend
    xlabel('Temperature (K)')
    ylabel('Lifetime (\mus)')
end




%% GN1881 contour plot
subplot(2,2,1)
%figure
contourLines = linspace(0,50,51);
contour(10.*GN1881Doping,-1000.*GN1881Defect,GN1881DChiSquare,contourLines)

hold on
plot(10.*[dopingInterests(1),dopingInterests(1)],-1000.*[defecthigh,defectlow],'k')
plot(10.*[dopingInterests(2),dopingInterests(2)],-1000.*[defecthigh,defectlow],'k')
plot(10.*[dopelow, dopehigh],-1000.*[defectInterests(1),defectInterests(1)],'k')
plot(10.*[dopelow, dopehigh],-1000.*[defectInterests(2),defectInterests(2)],'k')
plot(10.*[dopelow, dopehigh],-1000.*[defectInterests(3),defectInterests(3)],'k')


xlabel('Doping Density (x10^{14} cm^{-3})')
ylabel('Defect Level (meV)')
contourcbar
%% SAVE the doping Guesses, defect Guesses and DeltaChisquare to files
filePrefix = 'GN1881_C-S_';
txtString = strcat('DopingLow:','\t',num2str(dopelow),'\n',...
                   'DopingHigh:','\t',num2str(dopehigh),'\n',...
                   'DefectLow:','\t',num2str(defectlow),'\n',...
                   'DefectHigh:','\t', num2str(defecthigh),'\n',...
                   'Resolution doping:','\t',num2str(resDope),'\n',...
                   'Resolution defect:','\t',num2str(resDefect),'\n');              
fileID = fopen(strcat(filePrefix,'InputParameters.txt'),'w');
fprintf(fileID,txtString);
fclose(fileID);

writematrix(GN1881Defect, strcat(filePrefix,'DopingGuesses.csv'));
writematrix(GN1881Doping, strcat(filePrefix,'DefectGuesses.csv'));
writematrix(GN1881DChiSquare, strcat(filePrefix,'DeltaChiSquare.csv'));
save(strcat(filePrefix,'FunctionObjects'),'GN1881fobjects');

%% NEXT SAMPLE

%% Now do it for GN1880
filenameA = 'GN1880-combined10122020_LowestExcitation-12-11-2020_CarrierAnalysis';
filenameLtsA = 'GN1880-combined10122020_FittedResultsLT';
load(filenameA,'bestInputs','indexBest','fobjects')
load(filenameLtsA,'t0_vTemp','t0Sigma_vTemp','Temperatures')
% Doping in units of 1e15
dopelow = 0.03;
dopehigh = 1.5;
% Defect level in units of eV below the conduction band  which is the reference
defectlow = -0.105;
defecthigh = -0.075;
% Resolution of the fits
resDope = 70;
resDefect = 60;

[GN1880fobjects,GN1880Doping, GN1880Defect,GN1880DChiSquare] = lookAtChiSquare(filenameA,...
    filenameLtsA,...
    dopelow,dopehigh,defectlow,defecthigh,resDope, resDefect);


%% GN1880 contour plot
figure
contourLevelsGN1880 = linspace(0,50,51);
contour(GN1880Doping,-1000*GN1880Defect,GN1880DChiSquare, contourLevelsGN1880)
xlabel('Doping Concentration (x10^{15} cm^{-3})')
ylabel('Defect Level (meV)')
contourcbar
%% SAVE the doping Guesses, defect Guesses and DeltaChisquare to files
filePrefix = 'GN1880_';
txtString = strcat('DopingLow:','\t',num2str(dopelow),...
                   '\n','DopingHigh:','\t',num2str(dopehigh),...
                   '\n','DefectLow:','\t',num2str(defectlow),...
                   '\n','DefectHigh:','\t', num2str(defecthigh),...
                   '\n','Resolution doping:','\t',num2str(resDope),...
                   '\n','Resolution defect:','\t',num2str(resDefect));
               
fileID = fopen(strcat(filePrefix,'C-S_InputParameters.txt'),'w');
fprintf(fileID,txtString);
fclose(fileID);

writematrix(GN1880Defect, strcat(filePrefix,'C-S_DopingGuesses.csv'));
writematrix(GN1880Doping, strcat(filePrefix,'C-S_DefectGuesses.csv'));
writematrix(GN1880DChiSquare, strcat(filePrefix,'C-S_DeltaChiSquare.csv'));
save(strcat(filePrefix,'FunctionObjects'),'GN1880fobjects');

%% look for any negatives
[minimums, minIndices] = min(GN1881DChiSquare);
[mimimums2,minIndices2] = min(minimums);
minrow = minIndices(minIndices2);
figure
errorbar(Temperatures, t0_vTemp(1,:,end),t0Sigma_vTemp(1,:,end),'o','Displayname','Data')
hold on
plot(Temperatures,fobjects{indexBest}(Temperatures),'DisplayName','Best Fit (in paper)')
plot(Temperatures, GN1881fobjects{minrow,minIndices2}(Temperatures),'DisplayName','New Fit')
hold off
legend
xlabel('Temperature (K)')
ylabel('Lifetime (\mus)')
%% Fix doping but fit everythin else
load(filenameB);
load(filenameLtsB,'t0Sigma_vTemp');
bestFobject = fobjects{indexBest};
errors = bestInputs.Errors(2:3);

doping = 0.15;


% Calculate the best fit values from the software tool
yTheory = bestFobject(Temperatures);

% Calculate the best fit chi-square
chisquareOriginal  = sum( ((yTheory - t0_vTemp(1,:,end)')./t0Sigma_vTemp(1,:,end)').^2);
disp(strcat('Best fit Sample Chi square:', num2str(chisquareOriginal)));
% Now the reduced chi-square where it's the original chi sqaure divided by
% the degrees of freedom (N - parameters)
%reducedchiSquare = chisquareOriginal/(length(Temperatures) - length(coeffvalues(bestFobject)));

%prepare inputs for Matlab fit
xin = Temperatures';
yin = t0_vTemp(1,:,end)';
weights = 1./t0Sigma_vTemp(1,:,end).^2';
%weights = 1./t0Sigma_vTemp(1,:,end)'

%Change in chiSquare
DELTACHISQUARE = 0; 



problemParameters = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'p'];

problemValues = {bestFobject.a,bestFobject.b, bestFobject.c, bestFobject.d, bestFobject.e, bestFobject.f,...
                    bestFobject.g, bestFobject.h, bestFobject.k, bestFobject.l, bestFobject.m, bestFobject.n, doping};

ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
        'problem', problemParameters);

startPoint = [bestFobject.o,bestFobject.q, bestFobject.r];

% Run through all of the dopings
disp('Busy fitting');

       
fobjecta = fit(xin, yin, ft,...
    'problem', problemValues,...
    'Startpoint', startPoint,...
    'Weights', weights);
chisquarea = sum(((yin - fobjecta(xin)).^2.*weights));
DELTACHISQUAREa = chisquarea - chisquareOriginal;     
figure
disp('done!')
figure
 errorbar(Temperatures, t0_vTemp(1,:,end),t0Sigma_vTemp(1,:,end),'o','Displayname','Data')
    hold on
    plot(Temperatures,fobjects{indexBest}(Temperatures),'DisplayName','Best Fit (in paper)')
    plot(Temperatures, fobjecta(Temperatures),...
        'DisplayName','newFit')
xlabel('Temperature (K)')
ylabel('Lifetime (\mus)')
%legend()
%{
subplot(1,2,1)
 errorbar(Temperatures, t0_vTemp(1,:,end),t0Sigma_vTemp(1,:,end),'o','Displayname','Data')
    hold on
    plot(Temperatures,fobjects{indexBest}(Temperatures),'DisplayName','Best Fit (in paper)')
    plot(Temperatures, fobjectInterest{3}(Temperatures),...
        'DisplayName','newFit')

xlabel('Temperature (K)')
ylabel('Lifetime (\mus)')
%}



