function[fobjects, DOPE,DEFECT,DELTACHISQUARE] = lookAtChiSquare(carrierAnalysisFilename, lifetimeFilename, dopingLow,...
    dopingHigh, negativeDefectLow,negativeDefectHigh, resolutionDoping, resolutionDefect)

% Load the results extracted from the carrier analysis tool and explore the
% Chi-Sqaure metric as we vary the doping and the defect level

load(carrierAnalysisFilename);
load(lifetimeFilename,'t0Sigma_vTemp');
bestFobject = fobjects{indexBest};
errors = bestInputs.Errors(2:3);


dopingGuesses = linspace(dopingLow, dopingHigh, resolutionDoping);
defectGuesses = linspace(negativeDefectLow, negativeDefectHigh, resolutionDefect);



%Now create a mesh as this will end up in a 3d plot or countour
[DOPE, DEFECT] = meshgrid(dopingGuesses,defectGuesses);

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
%Output function handles
fobjects = cell(length(defectGuesses),length(dopingGuesses));
%Change in chiSquare
DELTACHISQUARE = NaN*ones(resolutionDefect, resolutionDoping);



problemParameters = ['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'p'; 'q'];

problemValues = {bestFobject.a,bestFobject.b, bestFobject.c, bestFobject.d, bestFobject.e, bestFobject.f,...
                    bestFobject.g, bestFobject.h, bestFobject.k, bestFobject.l, bestFobject.m, bestFobject.n, NaN, NaN};

ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r)',...
        'problem', problemParameters);

startPoint = [bestFobject.o, bestFobject.r];

% Run through all of the dopings
disp('Busy fitting');
for dopeIndex = 1:resolutionDoping
    % Run through all of the defect levels
    for defectIndex = 1:resolutionDefect
        %Initialize parameters
        p = DOPE(defectIndex, dopeIndex);
        q = DEFECT(defectIndex, dopeIndex);
        problemValues{13} = p;
        problemValues{14} = q;
       
        fobjects{defectIndex,dopeIndex} = fit(xin, yin, ft,...
            'problem', problemValues,...
            'Startpoint', startPoint,...
            'Weights', weights);
        chisquare = sum(((yin - fobjects{defectIndex,dopeIndex}(xin)).^2.*weights));
        DELTACHISQUARE(defectIndex,dopeIndex) = chisquare - chisquareOriginal;     
    end
end
disp('done!')


end