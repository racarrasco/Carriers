%% Load the carrier analysis results
%{
% weighting on the data
filenames = ["GN1886_LowestExcitation-Post-Rad-06-23-2021_CarrierAnalysis";...
    "GN1886_LowestExcitation-Post-Anneal-06-23-2021_CarrierAnalysis";...
    "C:\Users\rigoc\Desktop\KirtlandDocs\AFRL_DATA\paperData\InGaAs-InAsSbLifetimeManuscript(DONE)-2021\Submission3\Data\LT\GN1886\GN1886-combined10122020_LowestExcitation-03-17-2021-correctWeights_CarrierAnalysis"];
%}

% same weighting on pre-rad and post-rad
filenames = ["C:\Users\rigoc\Desktop\KirtlandDocs\AFRL_DATA\paperData\RadToleranceMinorityCarrierLifetime\DataAndResults\WeightingTempSweeps\GN1886_LowestExcitation-Post-Rad-06-23-2021_CarrierAnalysis";...
    "C:\Users\rigoc\Desktop\KirtlandDocs\AFRL_DATA\paperData\RadToleranceMinorityCarrierLifetime\DataAndResults\WeightingTempSweeps\GN1886_LowestExcitation-Post-Anneal-06-23-2021_CarrierAnalysis";...
    "C:\Users\rigoc\Desktop\KirtlandDocs\AFRL_DATA\paperData\InGaAs-InAsSbLifetimeManuscript(DONE)-2021\Submission3\Data\LT\GN1886\GN1886-combined10122020_LowestExcitation-03-17-2021-correctWeights_CarrierAnalysis"];



figureTitles = ["Post-Rad";"Post-Anneal"; "Pre-Rad"];
minimumBounds = ...
   [0, 10, 50;...
    0, 10, 30;...
    0, 8, 7];

%load(filename);

%{
 Read the function objects 
o = Auger overlap parameter |F1F2|
p = doping density in units of 1e15 cm^(-3)
q = defect level in reference to the conduction band
r = defect density in units of inverse meters
%}

%% Begin observing all of the best fits


numOfFits = length(filenames);

for fileind = 1:numOfFits
    load(filenames(fileind));
    
    figure('Name', figureTitles(fileind));
    %the number of fit results
    numOfFits = length(fobjects);


    % Instantiate the doping results
    dopings = nan.*ones(numOfFits,1); %In units of 1e15 cm^(-3)

    % Instantiate the rmse results
    rmses = nan.*ones(numOfFits,1);

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
        rmses(fitInd) = gofs(fitInd).rmse;

        %Take the defect levels
        defectLevels(fitInd) = fobjects{fitInd}.q;
    end
    subplot(2,2,1)
    plot(dopings,rmses, 'o')
    title(strcat(SampleName,' Best fits'))
    xlabel("doping, (x1e15 cm^{-3})")
    ylabel("RMSE")

    %plot the bloch overlap parameters
    subplot(2,2,2)
    plot(augers,rmses, 'o')
    xlabel("Bloch overlap")

    % plot the defect levels
    subplot(2,2,3)
    plot(defectLevels, rmses,'o')
    xlabel("Defect levels (meV)")
    ylabel("RMSE")
    
    
    subplot(2,2,4)
    plot(Temperatures, t0_vTemp(1,:,end),'o', 'DisplayName', 'Data');
    hold on
    for boundNumber = 1:3 
        %Take the rmses greater than the minimum values specified in
        %minimumBounds
        
        if(boundNumber ==3 && fileind ==3)
            %only 1 exception of the 9 other scenarios to look at
            dopingsGreaterIndices = dopings < minimumBounds(fileind, boundNumber);
        else
            dopingsGreaterIndices = dopings > minimumBounds(fileind, boundNumber);
            
            
            
        end
        
        
        
        rmsesGreater = rmses(dopingsGreaterIndices);
        fobjectsGreater = fobjects(dopingsGreaterIndices);
        [lowest, indLowest] = min(rmsesGreater);
        
        auger = fobjectsGreater{indLowest}.o;
        doping = fobjectsGreater{indLowest}.p;
        defectLevel = fobjectsGreater{indLowest}.q;
        defectDensity = fobjectsGreater{indLowest}.r;

        meStar = 0.026;
        mhStar = 0.3333;
        egTemp = varshniEinstein(Tprobe)/1000;
        eV = -egTemp;
        eC = 0;
        einf = 12.2;

        [nc, nv, ni, G] = calculateKeyParameters2(Tprobe, meStar, mhStar, egTemp,...
            einf);
        lt = calculateTotalLifetime(Tprobe, 'p', meStar, mhStar, egTemp, eV, eC,...
            einf, nc, nv, ni, G, phi, auger, doping, defectLevel, defectDensity);
        
        
      
       plot(Tprobe, lt, 'DisplayName', strcat('RMSE:', num2str(lowest,3)))
       xlabel("Temperature (K)")
       ylabel("Minority carrier lifetime (\mus)")
        
        
    end
    legend()
    hold off
    
        
    
    

end


