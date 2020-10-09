function[fobjects, gofs, outputs] =...
    fitLifetimeLink(xin, yin, lengths, type, eg, valenceEdge, conductionEdge,...
  einf, phi, meStar, mhStar, f1f2,...
  dopingDensity1,dopingDensity2, dopingDensity3,...
  defectLevel1,...
defectLevel3, defectDensity1, defectDensity2, defectDensity3, defectlevel2)
    % Calculate all of the lifetime components and add them together the rates
    % are added then the reciprocal of the total rate is the total lifetime
    % SEE THE calculateLifetimes.m file for meaning to the input parameters
    % perform a stochastic approach

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % yin   =  the data we will be performing on the fitting

    % These will be the fit parameters, and MATLAB arranges the fit parameters
    % alphabetically, not in the order in which they appear.

    % Give a weighting to each measurement and just use that value as the
    % weight
    % 1 us error is 1us, and 400 ns error is 400 ns etc.
    % w = yin';
    w = [];
    
    % initialize the Mersenne twister generator using a seed of 2
    s = rng(2);
    
    % Include the user guesses and size of the random guesses
    randomGuesses = 100; % This can change
    iterations = randomGuesses + 1; % Include the user guess
    guesses = zeros(iterations,  12);
    bestFits = zeros(iterations, 12);
    
    
    
    % Set the number of random guesses of the parameters and initialize the
    % array of random guesses
    guessInitialization = rand(randomGuesses, 12);
    
    %%% BEGIN STOCHASTIC GUESSING 
    %    first row is user input
    %    second row onward is stochastic guessing in all parameters
    %    columns are as follows: meStar, mhStar, f1f2, dopingDensity1,
    %    dopingDensity2, dopingDensity3, defectLevel1, defectLevel3,
    %    defectDensity1, defectDensity2, defectDensity3

    
    % meStar is between 0.01 and 1 m0
    meStarRange = [0.01, 1];
    guesses(2:end, 1) = meStarRange(1) + ...
        (meStarRange(end) - meStarRange(1)).*guessInitialization(:,1);

    %mhStar is between 0.2 and 110
    mhStarRange = [0.2, 110];
    guesses(2:end, 2) = mhStarRange(1) + ...
        (mhStarRange(end) - mhStarRange(1)).*guessInitialization(:,2);

    % f1f2 range between 0.1 and 0.3
    f1f2Range = [0.1, 0.3];
    guesses(2:end, 3) = f1f2Range(1)  + ...
        (f1f2Range(end) - f1f2Range(1)).*guessInitialization(:,3);

    % doping density between 1e13 and 1e17
    dopingDensityRange = [0.01, 100];
   
    guesses(2:end, 4:6) = dopingDensityRange(1) + ...
        (dopingDensityRange(end) - dopingDensityRange(1)).*guessInitialization(:,4:6);
    
    
    
     if(nargin ==20) 
    % [f1f2, dopingdensity1, dopingdensity2,
    % defectlevel, defectDensity1, defectdensity2, meStar, mhStar]
    guesses(1,1:end -1) = [meStar, mhStar, f1f2, dopingDensity1, dopingDensity2,...
        dopingDensity3, defectLevel1, defectLevel3, defectDensity1, defectDensity2,...
        defectDensity3];
    
    % defectLevel 1 is between 30 meV and 200 meV
    defectLevel1Range = [-.200, -.030];
    guesses(2:end, 7:8) = defectLevel1Range(1) + ...
        (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,7:8);

    % defect density 1 is between 0.1 m-1 and 60 m-1
    defectDensity1Range = [0.1, 60];
    guesses(2:end, 9:11) = defectDensity1Range(1) + ...
        (defectDensity1Range(end) - defectDensity1Range(1)).*guessInitialization(:,9:11);
     else
        
    % [f1f2, dopingdensity1, dopingdensity2,
    % defectlevel, defectDensity1, defectdensity2, meStar, mhStar]
    guesses(1,:) = [meStar, mhStar, f1f2, dopingDensity1, dopingDensity2,...
        dopingDensity3, defectLevel1,defectlevel2, defectLevel3,...
        defectDensity1, defectDensity2, defectDensity3];
    
    % defectLevel 1 is between 30 meV and 200 meV
    defectLevel1Range = [-.200, -.030];
    guesses(2:end, 7:9) = defectLevel1Range(1) + ...
        (defectLevel1Range(end) - defectLevel1Range(1)).*guessInitialization(:,7:9);

    % defect density 1 is between 0.1 m-1 and 60 m-1
    defectDensity1Range = [0.1, 60];
    guesses(2:end, 10:12) = defectDensity1Range(1) + ...
        (defectDensity1Range(end) - defectDensity1Range(1)).*guessInitialization(:,10:12); 
        
    end
   

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
    f1f2Limits = [0.01, 0.3];
    dopingLimits = [0.001, 1000]; %units of 1e15 cm-3
    defectLevelLimits = [-0.3,0]; %units of eV
    defectDensityLimits = [0,200]; 
    mhStarLimits = [0.2, 110];
    meStarLimits = [0.01, 1];
    
figure


for i = 1:iterations
    if(nargin==20)
    ft = fittype('calculateTotalLifetimeLink(x, a, b, c, d, e, f, g, h, k, l, m1, m2, m3, n1, n3, o1, o2, o3)',...
        'problem', ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g';]);
    startguess = guesses(i,1:end -1);
    
    [fobjects{i}, gofs(i), outputs(i)] = fit(xin, yin, ft,...
        'problem',{lengths, type, eg, valenceEdge, conductionEdge,einf, phi},...
        'Lower', [meStarLimits(1), mhStarLimits(1),f1f2Limits(1),...
        dopingLimits(1),dopingLimits(1),dopingLimits(1), defectLevelLimits(1),...
        defectLevelLimits(1), defectDensityLimits(1),defectDensityLimits(1),...
        defectDensityLimits(1)],...
        'Upper',[meStarLimits(2), mhStarLimits(2),f1f2Limits(2),...
        dopingLimits(2),dopingLimits(2),dopingLimits(2), defectLevelLimits(2),...
        defectLevelLimits(2), defectDensityLimits(2),defectDensityLimits(2),...
        defectDensityLimits(2)],...
        'Startpoint', startguess,...
        'Weights', w,...
        'MaxFunEvals', 2000,...
        'Maxiter', 200); 
    
    else
        ft = fittype('calculateTotalLifetimeLink2(x, a, b, c, d, e, f, g, h, k, l, m1, m2, m3, n1, n2, n3, o1, o2, o3)',...
        'problem', ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g';]);
    startguess = guesses(i,:);
    
    [fobjects{i}, gofs(i), outputs(i)] = fit(xin, yin, ft,...
        'problem',{lengths, type, eg, valenceEdge, conductionEdge,einf, phi},...
        'Lower', [meStarLimits(1), mhStarLimits(1),f1f2Limits(1),...
        dopingLimits(1),dopingLimits(1),dopingLimits(1), defectLevelLimits(1),...
        defectLevelLimits(1),defectLevelLimits(1), defectDensityLimits(1),defectDensityLimits(1),...
        defectDensityLimits(1)],...
        'Upper',[meStarLimits(2), mhStarLimits(2), f1f2Limits(2),...
        dopingLimits(2), dopingLimits(2), dopingLimits(2), defectLevelLimits(2),...
        defectLevelLimits(2), defectLevelLimits(2), defectDensityLimits(2),...
        defectDensityLimits(2), defectDensityLimits(2)],...
        'Startpoint', startguess,...
        'Weights', w,...
        'MaxFunEvals', 2000,...
        'MaxIter', 600); 
    end
    
    
    
    
    
    
    
    fobjecti = fobjects{i};
    yi = fobjecti(xin);
    
    plot(xin(1:lengths(1)), yin(1:lengths(1)),'o')
    hold on
    plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
        yin(lengths(1) + 1 :sum(lengths(1:2))),'+')
    lastIndex = sum(lengths(1:2));
    
    plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
        yin(lastIndex + 1: lastIndex+lengths(3)),'>');
    
    
    %Show the data that was shown
    plot(xin(1:lengths(1)), yi(1:lengths(1)),'k-')
    
    plot(xin(lengths(1) + 1:sum(lengths(1:2))),...
        yi(lengths(1) + 1 :sum(lengths(1:2))),'k--')
    lastIndex = sum(lengths(1:2));
    
    plot(xin(lastIndex + 1:lastIndex + lengths(3)),...
        yi(lastIndex + 1: lastIndex+lengths(3)),'k.');
    titleString = "Iteration: " + string(i) +"/" +...
        string(iterations)+ "; me/m0:  " + string(fobjecti.h) + ...
        "; mh/m0: " + string(fobjecti.k) + "; rmse: " + string(gofs(i).rmse);
    title(titleString);
    xlabel("Temperature (K)");
    ylabel("Minority carrier lifetime (us)");
   
    
    hold off
end
    
    

    
    
    
end