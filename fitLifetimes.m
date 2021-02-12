function[fobjects, gofs, outputs,inputParameters,bestIndex] = ...
 fitLifetimes(app, xin,yin, w, fitAugerOverlap,fixSRH1,fixSRHDefectLevel1, type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
  einf, Nc, Nv, ni,G, phi, f1f2, dopingDensity, defectLevel, defectDensity,...
  lowerN0, upperN0,lowerDefectLevel, upperDefectLevel,...
  lowerDefectDensity, upperDefectDensity, lowerF1F2, upperF1F2,...
  rngSeed, numberofGuesses,...
  defectLevel2, defectDensity2)

% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% SEE THE calculateLifetimes.m file for meaning to the input parameters

%perform a stochastic approach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%yin   =  the data we will be performing on the fitting

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, not in the order in which they appear.


% There are a lot of forks due to the multiple checkboxes
% So this is something that will be 


% initialize the Mersenne twister generator using a seed (default 2)
s = rng(rngSeed);
%Initialize rownames for best fit variables
fitNames = {};
if(nargin==33)
    defectLevel2 = 0;
    defectDensity2 = 0;
end
%Include the user guesses and size of the random guesses

randomGuesses = numberofGuesses; %this can change
iterations = randomGuesses + 1; %include the user guess
guesses = zeros(iterations,6); %guesses

% [f1f2, dopingdensity, defectlevel, defectDensity, defectlevel2,
% defectdensity2]
guesses(1,:) = [f1f2,dopingDensity, defectLevel, defectDensity,...
    defectLevel2, defectDensity2];


% Set the number of random guesses of the parameters and initialize the array 
% of random guesses

guessInitialization = rand(randomGuesses,6);

%%% BEGIN STOCHASTIC GUESSING 
%    first row is user input
%    second row onward is stochastic guessing in all parameters
%    columns are as follows: f1f2, dopingDensity, defectLevel1,
%    defectDensity1, defectLevel2, defectDensity2

% f1f2 range between 0.1 and 0.3
f1f2Range = [lowerF1F2*(1.1), upperF1F2*(0.90)];
guesses(2:end, 1) = f1f2Range(1)  +...
    (f1f2Range(end) - f1f2Range(1)).*guessInitialization(:,1);

% doping density between 7.25e13 and upperN0 - 0.45*upperN0
dopingDensityRange = [lowerN0*(1.1), upperN0*(0.90)];
guesses(2:end, 2) = dopingDensityRange(1) + ...
    (dopingDensityRange(end) - dopingDensityRange(1)).*guessInitialization(:,2);

% defectLevel 1 is between 30 meV and 300 meV
defectLevel1Range = [lowerDefectLevel*(0.90), upperDefectLevel*(0.9)];
guesses(2:end, 3) = defectLevel1Range(1) + ...
    (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,3);

% defect density 1 is between 0.1 m-1 and 60 m-1 (default)
defectDensity1Range = [lowerDefectDensity*(1.1), upperDefectDensity*(0.9)];
guesses(2:end, 4) = defectDensity1Range(1) + ...
    (defectDensity1Range(end) - defectDensity1Range(1)).*guessInitialization(:,4);

% defectLevel 2 has same range parameters as defectLevel 1 but different
% stochastic guesses
guesses(2:end, 5) = defectLevel1Range(1) + ...
    (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,5);


% defectDensity 2 has same range parameters, but guesses will be different
guesses(2:end, 6) = defectDensity1Range(1) + ...
    (defectDensity1Range(end) - defectDensity1Range(1)).*guessInitialization(:, 6);

%%% END STOCHASTIC GUESSING


%%% PREALLOCATE FITTING RESULTS
gofs  = repmat(struct('sse', 0,...
    'rsquare', 0,...
    'dfe',0,...
    'adjrsquare', 0,...
    'rmse',0), iterations, 1);


outputs = repmat(struct('numobs',0,'numparam', 0,...
    'residuals',[],...
    'Jacobian',0,...
    'exitflag',0,...
    'firstorderopt',0,...
    'iterations',0,...
    'funcCount',0,...
    'cgiterations',0,...
    'algorithm','algorithm',...
    'stepsize',0,...
    'message','message'), iterations, 1);
fobjects = cell(iterations,1);
%%% END PREALLOCATION

%%% BEGIN FITTING LIMITS SO FIT WILL NOT GO ABOVE OR BELOW THESE LIMITS
%%% SPECIFIED BY THE USER
f1f2Limits = [lowerF1F2, upperF1F2];
dopingLimits = [lowerN0, upperN0]; %units of 1e15 cm-3
defectLevelLimits = [lowerDefectLevel,upperDefectLevel]; %units of eV
defectDensityLimits = [lowerDefectDensity,upperDefectDensity];




casenumber = 0;
message = "";
maxFunevals = 1000;
maxIter = 600;



%%% INITIALIZE PROBLEM PARAMETERS
problemParameters = ['a'; 'b';'c';'d'; 'e'; 'f'; 'g';'h';'k';'l';'m';'n'];
%Initialized Problem Values
problemValues = {type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                    einf, Nc, Nv, ni, G, phi};
initialSize = size(problemParameters,1);

%%% ADD MORE PROBLEM PARAMETERS OR FIXED PARAMETERS
if(nargin ==33) %do we have an undoped sample or just 1 SRH level
    iterationGuesses = [guesses(:,1), guesses(:,2), guesses(:,3), guesses(:,4)];
    fitNames = {'f1f2','dopingDensity_1e15_cm3', 'defectLevel_eV','defectDensity_m1'};
    lower = [f1f2Limits(1), dopingLimits(1),defectLevelLimits(1),defectDensityLimits(1)];
    upper = [f1f2Limits(2), dopingLimits(2),defectLevelLimits(2),defectDensityLimits(2)];
    
    %Fix the SRH defect level
    if(fixSRHDefectLevel1) 
       %Add the problem parameter in 'problem' parameter
       problemParameters(initialSize + 1) = 'q';
       %Add the problem value in problemValue
       problemValues{initialSize + 1} = defectLevel;
       %increase the size of the problem array
       initialSize = initialSize + 1;

       %Remove SRH from array of guesses
       iterationGuesses(:,3) = [];
       fitNames(3) = [];
       lower(3) = [];
       upper(3) =[];

    end
    
    %fix the bloch overlap value
    if(~fitAugerOverlap) 
        %Add the block overlap parameter as a FIXED 'problem' parameter
        problemParameters(initialSize + 1) = 'o';
        %the user input is the problem vale
        problemValues{initialSize + 1} = f1f2;
        %expand the size
        initialSize = initialSize + 1;   
        %remove the bloch overlap from the guess and limits
        iterationGuesses(:,1) = [];
        fitNames(1) = [];
        lower(1) = [];
        upper(1) = [];
    end
    
    
    
    
    ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
        'problem', problemParameters);
else %We have a doped sample
    iterationGuesses = guesses;
    fitNames = {'f1f2','dopingDensity_1e15_cm3', 'defectLevel_eV','defectDensity_m1',...
        'defecLevel2_ev','defectDensity2_m1'};
    lower = [f1f2Limits(1), dopingLimits(1), defectLevelLimits(1), defectDensityLimits(1),...
        defectLevelLimits(1), defectDensityLimits(1)];
    upper = [f1f2Limits(2), dopingLimits(2), defectLevelLimits(2), defectDensityLimits(2),...
        defectLevelLimits(2), defectDensityLimits(2)];

    if(fixSRH1) %Fix the SRH values
        problemParameters(initialSize + 1) = 'q';
        problemValues{initialSize + 1} = defectLevel;
        problemParameters(initialSize + 2) = 'r';
        problemValues{initialSize + 2} = defectDensity;
        initialSize = initialSize + 2;

        iterationGuesses(:,4) = [];
        iterationGuesses(:,3) = [];

        fitNames(4) = [];
        fitNames(3) = [];

        lower(4) = [];
        lower(3) =[];

        upper(4) = [];
        upper(3) = [];

    end
    if(~fitAugerOverlap)% Fix the bloch overlap value
        problemParameters(initialSize + 1) = 'o';
        problemValues{initialSize + 1} = f1f2;
        initialSize = initialSize + 1;
        iterationGuesses(:,1) = [];
        fitNames(1) = [];
        lower(1) = [];
        upper(1) = [];
    end
    
   
        ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r,s,t)',...
        'problem', problemParameters);  
end






%%% BEGIN STOCHASTIC FITTING ALGORITHM
for i = 1: iterations
    [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                    'problem',problemValues,...
                    'Lower', lower,...
                    'Upper', upper, ...
                    'Startpoint', iterationGuesses(i,:),...
                    'Weights', w,...
                    'MaxFunEvals', maxFunevals,...
                    'MaxIter', maxIter);
    %PLOT the iterations of the fits on the app

    dopingDensitybest = fobjects{i}.p;
    defectLevelbest = fobjects{i}.q;
    defectDensitybest = fobjects{i}.r;
    blochOverlapbest = fobjects{i}.o;
    ColorWheel2 = [ "k--" ; "k-." ; "k:" ; "k-"; "-.r*"];
    

    if(nargin == 33)
    [totalLifetime, tauRad, tauSRH,tauAug] = calculateLifetimes(xin, type,meStar,...
                    mhStar, eg, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlapbest, dopingDensitybest, defectLevelbest, defectDensitybest);
    else
        defectLevel2best = fobjects{i}.s;
        defectDensity2best = fobjects{i}.t;
    [totalLifetime, tauRad, tauSRH, tauSRH2, tauAug] = calculateLifetimes2(xin, type,meStar,...
                    mhStar, eg, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlapbest, dopingDensitybest, defectLevelbest, defectDensitybest,...
                    defectLevel2best, defectDensity2best);
    end
        
                errorbar(app.UIAxes,xin, yin,1./w,'ko')
                hold(app.UIAxes, 'on')
                plot(app.UIAxes,xin,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                    
                plot(app.UIAxes,xin,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                if(nargin>33)
                    plot(app.UIAxes, xin, tauSRH2, ColorWheel2(5,:), 'Displayname','SRH2')
                    app.DefectLevel2EditField.Value = -defectLevel2best*1000;
                    app.DefectDensity2EditField.Value = defectDensity2best;
                end
                plot(app.UIAxes,xin,tauAug,ColorWheel2(3,:), 'Displayname', 'Auger')
                    plot(app.UIAxes,xin,totalLifetime, ColorWheel2(4,:),'Displayname','Total')
                    legend(app.UIAxes)
                    set(app.UIAxes, 'YScale', 'log')
                    
                    hold(app.UIAxes, 'off')
                    ylim(app.UIAxes,[1e-2, 1e3])
                    yticks(app.UIAxes,'auto')  
                    titleString = "Temperature Dependence (" + string(i) + "/" +...
                        string(iterations) + ")";
                    title(app.UIAxes,titleString);
                    app.DopingDensityEditField.Value = dopingDensitybest;
                    app.DefectLevel1EditField.Value = -1000*defectLevelbest;
                    app.DefectDensity1EditField.Value = defectDensitybest;
                    app.BlochoverlapEditField.Value = blochOverlapbest;  
    
end

% Get the best fit results from the guesses
[bestrmse,bestIndex] = min([gofs.rmse]);
%output results for the jacobian
bestOutput = outputs(bestIndex);
%best function object
bestfobject = fobjects{bestIndex};
%size of temperaure
tempsize = size(yin,2);
%errors
[sigmapBest,correlationBest] = calculateErrors(tempsize,bestOutput.Jacobian,...
    bestrmse,w);
coeffss = coeffvalues(bestfobject);
inputParameters = saveCoefficientsAndErrors(coeffss,fitNames, sigmapBest);          
      
end