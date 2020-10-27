function[fobjects, gofs, outputs,inputParameters,bestIndex] = ...
 fitLifetimes2(app, xin,yin,w, fitAugerOverlap,fixSRH1, type, eg, valenceEdge, conductionEdge,...
  einf, phi, meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity,lowerN0,upperN0,...
  lowerDefectLevel, upperDefectLevel,...
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

% Case 2: nargin == 28, but Bloch checkbox is false
% fit doping, defect1 level, and defect1 density but do not fit Bloch overlap 
%{
Don't need this anymore since the weights are an input now
% Give a weighting to each measurement and just use that value as the
% weight
% 1 us error is 1us,   and 400 ns error is 400 ns etc.
w = [];
%w = yin';
%}

% initialize the Mersenne twister generator using a seed of 2
s = rng(rngSeed);
%initialize rownames for best fit variables
fitNames = {};
if(nargin == 28)
   defectLevel2 = 0; 
   defectDensity2 = 0;

end
% Include the user guesses and size of the random guesses
randomGuesses = numberofGuesses; % this can change
iterations = randomGuesses + 1; % include the user guess
guesses = zeros(iterations,8);  % guesses
bestFits = zeros(iterations,8); % best fit parameters

% [f1f2, dopingdensity, defectlevel, defectDensity, defectlevel2,
% defectdensity2]
guesses(1,:) = [meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity,...
    defectLevel2, defectDensity2];

% Set the number of random guesses of the parameters and initialize the
% array of random guesses
guessInitialization = rand(randomGuesses, 8);

%%% BEGIN STOCHASTIC GUESSING INITIALIZATION 
%    first row is user input
%    second row onward is stochastic guessing in all parameters
%    columns are as follows: f1f2, dopingDensity, defectLevel1,
%    defectDensity1, defectLevel2, defectDensity2

% f1f2 range between 0.1 and 0.3
f1f2Range = [lowerF1F2*(1.1), upperF1F2*(0.90)];
guesses(2:end, 1) = f1f2Range(1)  +...
    (f1f2Range(end) - f1f2Range(1)).*guessInitialization(:,1);

% doping density between 1.45e14 and upperN0 - 0.45*upperN0
dopingDensityRange = [lowerN0*(1.1), upperN0*(0.90)];
guesses(2:end, 2) = dopingDensityRange(1) + ...
    (dopingDensityRange(end) - dopingDensityRange(1)).*guessInitialization(:,2);

% defectLevel 1 is between 30 meV and 300 meV
defectLevel1Range = [lowerDefectLevel*(0.90), upperDefectLevel(0.9)];
guesses(2:end, 3) = defectLevel1Range(1) + ...
    (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,3);

% defect density 1 is between 0.1 m-1 and 60 m-1
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

% meStar is between 0.01 and 1 m0
meStarRange = [0.01, 1];
guesses(2:end, 7) = meStarRange(1) + ...
    (meStarRange(end) - meStarRange(1)).*guessInitialization(:,7);

%mhStar is between 0.2 and 110
mhStarRange = [0.2, 110];
guesses(2:end, 8) = mhStarRange(1) + ...
    (mhStarRange(end) - mhStarRange(1)).*guessInitialization(:,8);

%%% END STOCHASTIC GUESSING INITIALIZATION


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
f1f2Limits = [0.01, 0.3];
dopingLimits = [0.05, upperN0]; %units of 1e15 cm-3
defectLevelLimits = [-0.3,0]; %units of eV
defectDensityLimits = [0,100];
mhStarLimits = [0.2, 110];
meStarLimits = [0.01, 1];

casenumber = 0;
message = "";
maxFunevals = 1000;
maxIter = 600;
%%% BEGIN STOCHASTIC FITTING ALGORITHM
%%% CONSIDER ALL POSSIBLE SWITCH CASES WITH FITTING EFFECTIVE MASSES
for i = 1:iterations
    if(nargin == 28) 
        if(fitAugerOverlap) 
                casenumber = 1;
                message = "Fit the effective masses, bloch overlap, doping, and SRH1 (6 parameters)";
                % CASE 1
                % Fit Effective masses auger overlap, doping,
                % defect level and defect density,
                ft = fittype('calculateTotalLifetime3(x, a, b, c, d, e, f, g, h, k, l, m, n)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f']);

                %guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity];
                guess = [guesses(i,7), guesses(i,8), guesses(i, 1),...
                    guesses(i, 2), guesses(i, 3) guesses(i, 4)];
                fitNames = {'meStar', 'mhStar', 'f1f2','dopingDensity_1e15_cm3',...
                    'defectLevel_eV','defectDensity_m1'};
                
                %Perform the fit with given parameters and lower/upper
                %limits, startpoint is given by "guess"
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, eg, valenceEdge, conductionEdge, einf, phi},...
                'Lower', [meStarLimits(1), mhStarLimits(1), f1f2Limits(1),...
                 dopingLimits(1), defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper', [meStarLimits(2), mhStarLimits(2), f1f2Limits(2),...
                 dopingLimits(2), defectLevelLimits(2), defectDensityLimits(2)],...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter',maxIter);
            
            bestFits(i,1:4) = [fobjects{i}.k, fobjects{i}.l, fobjects{i}.m,...
                fobjects{i}.n];
            bestFits(i,7:8) = [fobjects{i}.g, fobjects{i}.h];


        else
            casenumber = 2;
            message = "Fit effective masses, doping density, and SRH1 (5 parameters)";
            % CASE 2
            % Only fit 5 parameters which are effective masses (2), doping Density,
            % defect Level, and defect density,NOT bloch overlap
            ft = fittype('calculateTotalLifetime3(x, a, b, c, d, e, f, g, h, k, l, m, n)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k';]);
            % The guess is doping density, defect Level and defect density
            % guess = [meStar, mhStar,dopingDensity, defectLevel, defectDensity];
            guess = [guesses(i, 7), guesses(i, 8), guesses(i, 2), guesses(i, 3),...
                guesses(i, 4)];
            fitNames = {'meStar', 'mhStar', 'dopingDensity_1e15_cm3','defectLevel_eV',...
                'defectDensity_m1'};
            [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge,...
            einf, phi, f1f2},...
            'Lower', [meStarLimits(1), mhStarLimits(1), dopingLimits(1),...
            defectLevelLimits(1), defectDensityLimits(1)],...
            'Upper', [meStarLimits(2), mhStarLimits(2), dopingLimits(2),...
            defectLevelLimits(2), defectDensityLimits(2)],...
            'Startpoint', guess,...
            'Weights', w,...
            'MaxFunEvals', maxFunevals,...
            'MaxIter', maxIter);
        bestFits(i,2:4) = [fobjects{i}.l, fobjects{i}.m,...
                fobjects{i}.n];
        bestFits(i,7:8) = [fobjects{i}.g, fobjects{i}.h];
         

        end







    % We have more input arguments for fitting, which means that we have intentional
    % doping
    else
        if (fitAugerOverlap)
            if(fixSRH1) 
            casenumber = 3;
            message = "Fit effective masses, bloch overlap, doping, and SRH2";
            %FIX SRH1 parameters, fit effective masses
            %SRH2 parameters, Bloch overlap parameter,
            %doping (Fit 6 parameters)
                ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'm'; 'n']);

            %Fit parameter effective masses (g,h) k(Bloch overlap),
            %l (doping density), o(defectLevel2), p(defectDensity2) 

                %guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel2, defectDensity2];
                
                guess = [guesses(i,7), guesses(i,8), guesses(i,1), guesses(i,2),...
                    guesses(i,5), guesses(i,6)];
                fitNames = {'meStar','mhStar', 'f1f2','dopingDenstiy_1e15_cm3',...
                    'defectLevel2_eV','defectDenstiy2_m1'};


                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, eg, valenceEdge, conductionEdge,...
                einf, phi,defectLevel, defectDensity},...
                'Lower', [meStarLimits(1), mhStarLimits(1), f1f2Limits(1),...
                dopingLimits(1), defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper', [meStarLimits(2), mhStarLimits(2), f1f2Limits(2),...
                dopingLimits(2), defectLevelLimits(2), defectDensityLimits(2)], ...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals',maxFunevals,...
                'MaxIter',maxIter); 
            bestFits(i,1:2) = [fobjects{i}.k, fobjects{i}.l];
            bestFits(i,5:end) = [fobjects{i}.o, fobjects{i}.p, fobjects{i}.g,...
                fobjects{i}.h];
            else
            casenumber = 4;
            message = "Fit effective masses, bloch overlap, doping density, SRH1, and SRH2";
            % Fit SRH1, SRH2, Bloch overlap parameter, doping and effective maases (8 parameters)
                ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f']);

            %Fit parameter effective mass (g,h) k(Bloch overlap) l(doping density)
            % m(defectLevel1) n(defectDensity1) o(defect level2) p(defect density2)
              %  guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity,...
              %     defectLevel2, defectDensity2];
              
                guess = [guesses(i,7), guesses(i,8), guesses(i,1), guesses(i,2),...
                    guesses(i,3), guesses(i, 4), guesses(i,5), guesses(i,6)];
                fitNames = {'meStar','mhStar','f1f2','dopingDensity_1e15_cm3',...
                    'defectLevel_eV','defectDensity_m1','defectLevel2_eV','defectDensity2_m1'};
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, eg, valenceEdge, conductionEdge,...
                einf, phi},...
                'Lower', [meStarLimits(1), mhStarLimits(1), f1f2Limits(1),...
                dopingLimits(1), defectLevelLimits(1), defectDensityLimits(1),...
                defectLevelLimits(1), defectDensityLimits(1)],...
                ...
                'Upper', [meStarLimits(2), mhStarLimits(2), f1f2Limits(2),...
                dopingLimits(2), defectLevelLimits(2), defectDensityLimits(2),...
                defectLevelLimits(2), defectDensityLimits(2)],...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals',maxFunevals,...
                'MaxIter',maxIter);
            bestFits(i,:) =[fobjects{i}.k, fobjects{i}.l, fobjects{i}.m,...
                fobjects{i}.n, fobjects{i}.o, fobjects{i}.p,...
                fobjects{i}.g, fobjects{i}.h];



            end


        else 
            if(fixSRH1)
                casenumber = 5;
                message = "Fit effective masses, doping and SRH2";
                %FIX SRH1 and Bloch overlap. 
                %FIX SRH1 parameters and Bloch overlap, fit SRH2 parameters,
                %doping and effective masses (Fit 5 parameters)
                ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k'; 'm'; 'n']);

            %Fit parameter  p (doping density), s(defectLevel2), t(defectDensity2) 
                %guess = [meStar, mhStar, dopingDensity, defectLevel2, defectDensity2];
                guess = [guesses(i,7), guesses(i,8), guesses(i,2), guesses(i,5),...
                    guesses(i,6)];
                
                fitNames = {'meStar','mhStar','dopingDensity_1e15_cm3','defectLevel2_eV',...
                    'defectDensity_m1'};

                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, eg, valenceEdge, conductionEdge,...
                einf,  phi,f1f2, defectLevel, defectDensity},...
                'Lower', [meStarLimits(1), mhStarLimits(1), dopingLimits(1),...
                defectLevelLimits(1), defectDensityLimits(1)],...
                'Upper',[meStarLimits(2), mhStarLimits(2), dopingLimits(2),...
                defectLevelLimits(2), defectDensityLimits(2)],...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter', maxIter);
            bestFits(i,2) = [fobjects{i}.l];
            bestFits(i,5:6) = [fobjects{i}.o, fobjects{i}.p];
            bestFits(i,7:8) = [fobjects{i}.g, fobjects{i}.h];
            else
                casenumber = 6;
                message = "Fit effective masses, doping, SRH1 and SRH2";
                % Fit effective mass, SRH1, SRH2, doping, FIX Bloxh overlap parameter (7 parameters)
                ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
                'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k']);
            
            %Fit parameter gh (effective masses) l(doping density) m(defectLevel1) n(defectDensity1)
            % o(defect level2) (defect density2)
                guess = [guesses(i,7), guesses(i,8), guesses(i,2), guesses(i,3), guesses(i,4),...
                    guesses(i,5), guesses(i,6)];
                fitNames = {'meStar', 'mhStar','dopingDensity_1e15_cm3','defectLevel_eV','defectDensity_m1',...
                    'defectLevel2_eV','defectDensity2_m1'};
                [fobjects{i}, gofs(i), outputs(i)] = fit(xin',yin',ft,...
                'problem',{type, eg, valenceEdge, conductionEdge,...
                einf, phi, f1f2},...
                'Lower', [meStarLimits(1), mhStarLimits(1), dopingLimits(1),...
                 defectLevelLimits(1), defectDensityLimits(1),...
                 defectLevelLimits(1), defectLevelLimits(1)],...
                'Upper', [meStarLimits(2), mhStarLimits(2), dopingLimits(2),...
                 defectLevelLimits(2), defectDensityLimits(2),...
                 defectLevelLimits(2), defectLevelLimits(2)],...
                'Startpoint', guess,...
                'Weights', w,...
                'MaxFunEvals', maxFunevals,...
                'MaxIter',maxIter);
            bestFits(i, 2:8) = [fobjects{i}.l, fobjects{i}.m, fobjects{i}.n,...
                fobjects{i}.o, fobjects{i}.p, fobjects{i}.g, fobjects{i}.h];
            
            end



        end
    end
     
    %PLOT the iterations of the fits on the app

    dopingDensitybest = fobjects{i}.l;
    defectLevelbest = fobjects{i}.m;
    defectDensitybest = fobjects{i}.n;
    blochOverlapbest = fobjects{i}.k;
    meStarBest = fobjects{i}.g;
    mhStarBest = fobjects{i}.h;
    ColorWheel2 = ["k--"; "k-."; "k:"; "k-"; "-.r*"];
    [Nc, Nv, ni, G] = calculateKeyParameters2(xin, meStarBest, mhStarBest,...
    eg, einf);

    if(nargin == 28)
    [totalLifetime, tauRad, tauSRH,tauAug] = calculateLifetimes(xin, type,meStarBest,...
                    mhStarBest, eg, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlapbest, dopingDensitybest, defectLevelbest, defectDensitybest);
    else
        defectLevel2best = fobjects{i}.o;
        defectDensity2best = fobjects{i}.p;
    [totalLifetime, tauRad, tauSRH, tauSRH2, tauAug] = calculateLifetimes2(xin, type,meStarBest,...
                    mhStarBest, eg, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlapbest, dopingDensitybest, defectLevelbest, defectDensitybest,...
                    defectLevel2best, defectDensity2best);
    end
        
                errorbar(app.UIAxes,xin, yin,1./w,'ko')
                hold(app.UIAxes, 'on')
                plot(app.UIAxes,xin,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                    
                plot(app.UIAxes,xin,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                app.DopingDensityEditField.Value = dopingDensitybest;
                app.DefectLevel1EditField.Value = -defectLevelbest*1000;
                app.DefectDensity1EditField.Value = defectDensitybest;
                if(nargin>28)
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
                app.mhEditField.Value = mhStarBest;
                app.meEditField.Value = meStarBest;
                app.BlochoverlapEditField.Value = blochOverlapbest;        
end
%Get the best fit results from the 101 guesses
[bestrmse, bestIndex] = min([gofs.rmse]);

%output resutls for the jacobian
bestOutput = outputs(bestIndex);
%best function object
bestfobject = fobjects{bestIndex};
%size of temperature
tempsize = size(yin,2);
%errors
[sigmapBest, correlationBest] = calculateErrors(tempsize, bestOutput.Jacobian,...
    bestrmse,w);
coeffss = coeffvalues(bestfobject);
inputParameters = saveCoefficientsAndErrors(coeffss, fitNames,sigmapBest);




end