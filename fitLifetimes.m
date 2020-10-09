function[fobjects, gofs, outputs,inputParameters,bestIndex] = ...
 fitLifetimes(app, xin,yin, fitAugerOverlap,fixSRH1, type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
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
% Case 1: fit Bloch overlap, doping, defect1 level and defect1 density,
% there is no secondary SRH component

% Case 2: nargin == 31, but Bloch checkbox is false
% fit doping, defect1 level, and defect1 density but do not fit Bloch overlap 

% Give a weighting to each measurement and just use that value as the
% weight
% 1 us error is 1us,   and 400 ns error is 400 ns etc.
%w = yin';
w = [];

% initialize the Mersenne twister generator using a seed (default 2)
s = rng(rngSeed);
%Initialize rownames for best fit variables
fitNames = {};
if(nargin==31)
    defectLevel2 = 0;
    defectDensity2 = 0;
end
%Include the user guesses and size of the random guesses

randomGuesses = numberofGuesses; %this can change
iterations = randomGuesses + 1; %include the user guess
guesses = zeros(iterations,6); %guesses
bestFits = zeros(iterations,6);% best fit parameters

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
defectLevel1Range = [lowerDefectLevel*(90), upperDefectLevel*(1.1)];
guesses(2:end, 3) = defectLevel1Range(1) + ...
    (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,3);

% defect density 1 is between 0.1 m-1 and 60 m-1 (default)
defectDensity1Range = [lowerDefectDensity*1.1, upperDefectDensity*2];
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

%%% BEGIN FITTING LIMITS
f1f2Limits = [lowerF1F2, upperF1F2];
dopingLimits = [lowerN0, upperN0]; %units of 1e15 cm-3
defectLevelLimits = [lowerDefectLevel,0]; %units of eV
defectDensityLimits = [lowerDefectDensity,upperDefectDensity];




casenumber = 0;
message = "";
maxFunevals = 1000;
maxIter = 600;

%%% BEGIN STOCHASTIC FITTING ALGORITHM
%%% CONSIDER ALL POSSIBLE SWITCH CASES WITH FIXED EFFECTIVE MASS
for i = 1: iterations
    if(nargin == 31)
        if(fitAugerOverlap)
                casenumber = 1;
                message = "Fit the Bloch overlap, doping, defectLevel1 and defect density product (4 parameters)";
                % CASE 1
                % Fit auger overlap, doping, defect level and defect density
                ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n']);

                %guess = [f1f2, dopingDensity, defectLevel, defectDensity]; 
                guess = [guesses(i,1), guesses(i,2), guesses(i,3), guesses(i,4)];
                fitNames = {'f1f2','dopingDensity_1e15_cm3', 'defectLevel_eV','defectDensity_m1'};


                 [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                    'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                    einf, Nc, Nv, ni, G, phi},...
                    'Lower', [f1f2Limits(1) , dopingLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1)],...
                    'Upper', [f1f2Limits(2) , dopingLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2)], ...
                    'Startpoint', guess,...
                    'Weights', w,...
                    'MaxFunEvals', maxFunevals,...
                    'MaxIter', maxIter);  
                bestFits(i,1:4) = [fobjects{i}.o, fobjects{i}.p, fobjects{i}.q,...
                    fobjects{i}.r];
        else
            casenumber = 2;
            message = "Fit the doping, defect Level  and defect density product (3 parameters)";
            % CASE 2
            % Only fit 3 parameters which are doping Density, defect Level, and
            % defect density, NOT bloch overlap

            ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o']);
            % The guess is doping density, defect Level and defect density
            %guess = [dopingDensity, defectLevel, defectDensity];
            guess = [guesses(i,2), guesses(i, 3), guesses(i,4)];

            %Names of the fitting parameters will come in handy when saving
            %the best iteration
            fitNames = {'dopingDensity_1e15_cm3', 'defectLevel_eV','defectDensity_m1'};


            [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi, f1f2},...
            'Lower', [dopingLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1)],...
            'Upper', [dopingLimits(2),defectLevelLimits(2),...
            defectDensityLimits(2)], ...
            'Startpoint', guess,...
            'Weights', w,...
            'MaxFunEvals', maxFunevals,...
            'MaxIter', maxIter);
            bestFits(i,2:4) = [fobjects{i}.p, fobjects{i}.q,...
                    fobjects{i}.r];

        end 


    % We have more input arguments for fitting, which means that we have intentional
    % doping
    else
        if (fitAugerOverlap)
            if(fixSRH1) 
            casenumber = 3;
            message = "Fit the Bloch overlap, doping, defectLevel1 and defect density product (4 parameters)";
            % CAES 3
            %FIX SRH1 parameters, fit SRH2 parameters, Bloch overlap parameter,
            %doping (Fit 4 parameters)
                ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'q'; 'r']);

            %Fit parameter o(Bloch overlap), p (doping density), s(defectLevel2), t(defectDensity2) 
                % guess = [f1f2, dopingDensity, defectLevel2, defectDensity2];
                 guess = [guesses(i,1), guesses(i,2), guesses(i,5), guesses(i,6)];
                 
                fitNames = {'f1f2','dopingDensity_1e15_cm3', 'defectLevel2_eV','defectDensity2_m1'};

                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                einf, Nc, Nv, ni, G, phi,defectLevel, defectDensity},...
                'Lower', [f1f2Limits(1) , dopingLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper', [f1f2Limits(2) , dopingLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2)], ...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter', maxIter);
                bestFits(i,1:2) = [fobjects{i}.o, fobjects{i}.p];
                bestFits(i,5:6) = [fobjects{i}.s, fobjects{i}.t];
            
            else
            casenumber = 4;
            message = "Fit the Bloch overlap, doping, SRH1 and SRH2 (6 parameters)";
            % Fit SRH1, SRH2, Bloch overlap parameter, doping (6 parameters)
                ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n']);

            %Fit parameter o(Bloch overlap) p(doping density) q(defectLevel1) r(defectDensity1)
            % s(defect level2) t(defect density2)
                %guess = [f1f2, dopingDensity, defectLevel, defectDensity,...
                    %defectLevel2, defectDensity2];
                guess = guesses(i,:);
                fitNames = {'f1f2','dopingDensity_1e15_cm3', 'defectLevel_eV','defectDensity_m1',...
                    'defectLevel2_eV','defectDensity2_m1'};
                
                
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                einf, Nc, Nv, ni, G, phi},...
                'Lower', [f1f2Limits(1), dopingLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper', [f1f2Limits(2),  dopingLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2)],...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter', maxIter);
                bestFits(i,1:end) = [fobjects{i}.o, fobjects{i}.p, fobjects{i}.q,...
                    fobjects{i}.r, fobjects{i}.s, fobjects{i}.t];
            end


        else 
            if(fixSRH1)
                casenumber = 5;
                message = "Fit the doping, SRH2 (3 parameters)";
                 %FIX SRH1 and Bloch overlap. 
                 %FIX SRH1 parameters and Bloch overlap, fit SRH2 parameters,
                 %doping (Fit 3 parameters)
                ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o'; 'q'; 'r']);

                
            %Fit parameter  p (doping density), s(defectLevel2), t(defectDensity2) 
                 %guess = [dopingDensity, defectLevel2, defectDensity2];
                guess = [guesses(i,2), guesses(i,5), guesses(i,6)];
                fitNames = {'dopingDensity_1e15_cm3','defectLevel2_eV','defectDensity2_m1'};
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                einf, Nc, Nv, ni, G, phi,f1f2, defectLevel, defectDensity},...
                'Lower', [dopingLimits(1),defectLevelLimits(1),...
                defectDensityLimits(1)],...
                'Upper', [dopingLimits(2),defectLevelLimits(2),...
                defectDensityLimits(2)], ...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter', maxIter);
                bestFits(i,2) = [fobjects{i}.p];
                bestFits(i,5:6) = [fobjects{i}.s, fobjects{i}.t];
            else
                casenumber = 6;
                message = "Fit the SRH1, SRH2, and doping, (5 parameters)";
                % Fit SRH1, SRH2, doping, FIX Bloxh overlap parameter (5 parameters)
                ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n';'o']);
            % Fit parameter  p(doping density) q(defectLevel1) r(defectDensity1)
            % s(defect level2) t(defect density2)
                %guess = [ dopingDensity, defectLevel, defectDensity,...
                    %defectLevel2, defectDensity2];
                guess = [guesses(i, 2), guesses(i, 3) guesses(i,4), guesses(i,5)...
                    guesses(i,6)];
                
                fitNames = {'dopingDensity_1e15_cm3','defectLevel_eV','defectDensity_m1',...
                    'defectLevel2_eV','defectDensity2_m1'};
                
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
                einf, Nc, Nv, ni, G, phi, f1f2},...
                'Lower', [dopingLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1),...
                    defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper', [dopingLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2),...
                    defectLevelLimits(2), defectDensityLimits(2)], ...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter', maxIter);
                bestFits(i,2:end) = [fobjects{i}.p, fobjects{i}.q,...
                    fobjects{i}.r, fobjects{i}.s, fobjects{i}.t];
            end



        end

    end
    
    %PLOT the iterations of the fits on the app

    dopingDensitybest = fobjects{i}.p;
    defectLevelbest = fobjects{i}.q;
    defectDensitybest = fobjects{i}.r;
    blochOverlapbest = fobjects{i}.o;
    ColorWheel2 = [ "k--" ; "k-." ; "k:" ; "k-"; "-.r*"];
    

    if(nargin == 31)
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
        
                plot(app.UIAxes,xin, yin,'ko')
                hold(app.UIAxes, 'on')
                plot(app.UIAxes,xin,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                    
                plot(app.UIAxes,xin,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                if(nargin>31)
                    plot(app.UIAxes, xin, tauSRH2, ColorWheel2(5,:), 'Displayname','SRH2')
                    app.DefectLevel2EditField.Value = -defectLevel2best*1000;
                    app.DefectDensity2EditField.Value = defectDensity2best;
                end
                plot(app.UIAxes,xin,tauAug,ColorWheel2(3,:), 'Displayname', 'Auger')
                    plot(app.UIAxes,xin,totalLifetime, ColorWheel2(4,:),'Displayname','Total')
                    legend(app.UIAxes)
                    set(app.UIAxes, 'YScale', 'log')
                    
                    hold(app.UIAxes, 'off')
                    ylim(app.UIAxes,[1e-1, 1e3])
                    titleString = "Temperature Dependence (" + string(i) + "/" +...
                        string(iterations) + ")";
                    title(app.UIAxes,titleString);
                    app.DopingDensityEditField.Value = dopingDensitybest;
                    app.DefectLevel1EditField.Value = -1000*defectLevelbest;
                    app.DefectDensity1EditField.Value = defectDensitybest;
                    app.BlochoverlapEditField.Value = blochOverlapbest;
   
        
    
    
    
    
    
end
%%% What parameters were being fit?
indices = logical([0, 0, 0, 0, 0, 0]);
columns = {'F1F2', 'dopingDensity_1e15_cm3', 'DefectLevel_meV', 'defectDensity_m1',...
    'DefectLevel2_meV','DefectDensity2_m1'};
    switch(casenumber)
        case 1 %bloch overlap, doping, defectlevel1 and defectdensity
            indices = logical([1, 1, 1, 1, 0, 0]);   
        case 2
            indices = logical([0, 1, 1, 1, 0, 0]);
        case 3
            indices = logical([1, 1, 0, 0, 1, 1]);
        case 4
            indices = logical([1, 1, 1, 1, 1, 1]);
        case 5
            indices = logical([0, 1, 0, 0, 1, 1]);
        case 6
            indices = logical([0, 1, 1, 1, 1, 1]);
            
    end
% Get the best fit results from the 101 guesses
[bestrmse,bestIndex] = min([gofs.rmse]);
%output results for the jacobian
bestOutput = outputs(bestIndex);
%best function object
bestfobject = fobjects{bestIndex};
%size of temperaure
tempsize = size(yin,2);
%errors
[sigmapBest,correlationBest] = calculateErrors(tempsize,bestOutput.Jacobian,...
    bestrmse);
coeffss = coeffvalues(bestfobject);
inputParameters = saveCoefficientsAndErrors(coeffss,fitNames, sigmapBest);

            
            
      
end