%% Carrier Analysis viewer
% Load a ####_CarrierAnalysis.mat file and read the contents

filename = "GN2120-Combined2120LTvTemp2022_CarrierAnalysis";
load(filename);

%{
 Read the function objects 
o = Auger overlap parameter |F1F2|
p = doping density in units of 1e15 cm^(-3)
q = defect level in reference to the conduction band
r = defect density in units of inverse meters
%}

%% Begin observing all of the best fits

%the number of fit results
numOfFits = length(fobjects);

% Instantiate the doping results
dopings = nan.*ones(numOfFits,1); %In units of 1e15 cm^(-3)

% Instantiate the rmse results
rmsesAll = nan.*ones(numOfFits,1);

% Instantiate the defect level results
defectLevels = nan.*ones(numOfFits,1); %In units of meV

% Instantiate the Auger overlap
augers = nan.*ones(numOfFits,1); %In units of meV

for fitInd = 1:numOfFits
   %Take the Bloch overlap
   augers(fitInd)= fobjects{fitInd}.o;

   %Take the dopings
   dopings(fitInd) = fobjects{fitInd}.p;
   
   %Take the Root-mean-square
   rmsesAll(fitInd) = gofs(fitInd).rmse;
   
   %Take the defect levels
   defectLevels(fitInd) = fobjects{fitInd}.q;
end
subplot(2,2,1)
plot(dopings,rmsesAll, 'o')
title(strcat(SampleName,' Best fits'))
xlabel("doping, (x1e15 cm^{-3})")
ylabel("RMSE")

%plot the bloch overlap parameters
subplot(2,2,2)
plot(augers,rmsesAll, 'o')
xlabel("Bloch overlap")

% plot the defect levels
subplot(2,2,3)
plot(defectLevels, rmsesAll,'o')
xlabel("Defect levels (meV)")
ylabel("RMSE")

%% Find the global/local minimums in the best fit values


%Choose lower doping limit 0 will provide the same result as the Carrier
%analysis value

%Chossing a doping value just above the value with the lowest RMSE, then
%the following code will provide the next best result. To look at all of
%the minimums, just continue to increase this value above the subsequent
%doping values... 
% EX: with "dopinglowerLimit" =  0, you can have 4 local minimums at 1e15,
% 3e15 and 7e15. and then in return you will have the parameter results for
% the 1e15..If you set the "dopinglowerLimit" = 2, then you will have the
% best fit parameters at the 3e15 value

dopinglowerLimit = 0; %units of 1e15 cm^(-3) 


dopes = nan*ones(numOfFits,1);
rmses = nan*ones(numOfFits,1);

%read through dopings and pick values above dopinglowerLimit
for index = 1:numOfFits
    dopingValue = fobjects{index}.p;   
    if(dopingValue > dopinglowerLimit)     
        dopes(index)  = dopingValue;    
        rmses(index) = gofs(index).rmse;
    end   

end

%plot doping values above dopinglowerLimit
figure()
plot(dopes, rmses,'o')
xlabel('Doping density (x 10^{15}cm^{-3})')
ylabel('RMSE')
title(SampleName)

%Take the result with the minimum rmse above lowerdopinglimit
[value, index] = min(rmses);

%Extract parameters and prepare plot
fobjectinterest = fobjects{index};
dopetype = fobjectinterest.a;
mestr = fobjectinterest.b;
mhstr = fobjectinterest.c;
egtemp = varshniEinstein(Tprobe)/1000;
vedge = -egtemp;
cedge = 0;
einff = fobjectinterest.g;
nc = densityInConduction(Tprobe, mestr);
nv = densityInValence(Tprobe,mhstr);
ni = calculateIntrinsic(nc, nv, egtemp,Tprobe);
gAnalyt = calcAnalyticG(Tprobe,mestr, mhstr,ni, egtemp, einff);
phi =  fobjectinterest.n;
f1f2 = fobjectinterest.o;
dope = fobjectinterest.p;
defectlevel = fobjectinterest.q;
defectdens = fobjectinterest.r;


[total, rad, srh, auger] = calculateLifetimes(Tprobe,dopetype, mestr, mhstr,...
    egtemp, vedge, cedge, einff, nc, nv, ni, gAnalyt, phi, f1f2, dope,...
    defectlevel, defectdens);
figure
hold on
plot(Tprobe, srh,'-.','DisplayName','SRH')
plot(Tprobe, rad, '--','DisplayName','Radiative')
plot(Tprobe, auger,'.','DisplayName','Auger')
plot(Tprobe, total,'-','DisplayName','Total')
errorbar(Temperatures, t0_vTemp(1,:,end),t0Sigma_vTemp(1,:,end),'o',...
    'DisplayName','Data')
xlabel('Temperature (K)')
ylabel('Lifetime (\mus)')
titleName = sprintf('N_D = %0.2f \\times 10^{15} cm^{-3}; (E_C-E_t) =  %0.1f  meV; |F_1F_2| = %0.4f; \\sigmaN_t = %0.2f m^{-1}',...
    dope, defectlevel*1000, f1f2, defectdens);
title(titleName);

legend
