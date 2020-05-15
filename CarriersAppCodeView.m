classdef Carriers < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        DataFittingStatusEditFieldLabel  matlab.ui.control.Label
        ToolStatus                      matlab.ui.control.EditField
        UIAxes                          matlab.ui.control.UIAxes
        InputparametersandusercontrolPanel  matlab.ui.container.Panel
        MaterialDropDownLabel           matlab.ui.control.Label
        MaterialDropDown                matlab.ui.control.DropDown
        IntentionaldopingCheckBox       matlab.ui.control.CheckBox
        DopingTypeDropDownLabel         matlab.ui.control.Label
        DopingDropDown                  matlab.ui.control.DropDown
        DopingDensityEditFieldLabel     matlab.ui.control.Label
        DopingDensityEditField          matlab.ui.control.NumericEditField
        X1e15cm3Label                   matlab.ui.control.Label
        DefectLevel1EditFieldLabel      matlab.ui.control.Label
        DefectLevel1EditField           matlab.ui.control.NumericEditField
        meVLabel_7                      matlab.ui.control.Label
        FixCheckBox                     matlab.ui.control.CheckBox
        DefectDensity1EditFieldLabel    matlab.ui.control.Label
        DefectDensity1EditField         matlab.ui.control.NumericEditField
        m1Label_7                       matlab.ui.control.Label
        DefectLevel2EditFieldLabel      matlab.ui.control.Label
        DefectLevel2EditField           matlab.ui.control.NumericEditField
        meVLabel_8                      matlab.ui.control.Label
        DefectDensity2EditFieldLabel    matlab.ui.control.Label
        DefectDensity2EditField         matlab.ui.control.NumericEditField
        m1Label_8                       matlab.ui.control.Label
        BlochoverlapEditFieldLabel      matlab.ui.control.Label
        BlochoverlapEditField           matlab.ui.control.NumericEditField
        FitCheckBox                     matlab.ui.control.CheckBox
        ActiveRegionThicknessEditFieldLabel  matlab.ui.control.Label
        ActiveRegionThicknessEditField  matlab.ui.control.NumericEditField
        umLabel                         matlab.ui.control.Label
        CaplayerrefractiveindexatbandgapenergyEditFieldLabel  matlab.ui.control.Label
        CaplayerrefractiveindexatbandgapenergyEditField  matlab.ui.control.NumericEditField
        DataFittingPanel                matlab.ui.container.Panel
        TRPLFittedresultsMATfileLabel   matlab.ui.control.Label
        TRPLMATfilePath                 matlab.ui.control.EditField
        PLFittedresultsMATfileEditFieldLabel  matlab.ui.control.Label
        PLMATfilePath                   matlab.ui.control.EditField
        StartButton                     matlab.ui.control.Button
        SimulateButton                  matlab.ui.control.Button
        FitButton                       matlab.ui.control.Button
        SaveButton                      matlab.ui.control.Button
        AbortButton                     matlab.ui.control.Button
        NoPLDataCheckBox                matlab.ui.control.CheckBox
        OutputfilesuffixLabel           matlab.ui.control.Label
        UserOutputFileSuffix            matlab.ui.control.EditField
        ShowResultsCheckBox             matlab.ui.control.CheckBox
        ErrorEditFieldLabel             matlab.ui.control.Label
        ErrorReport                     matlab.ui.control.EditField
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: StartButton
        function Initialize(app, event)
            
            % Reset the error report and background return
            app.ErrorReport.Value = '';
            
            % Empty vectors from workspace that could confuse the tool if subsequent load steps fail
            Carrier_WorkspaceReset
            % Global variables for Previous/Repeat/Advance/Abort button presses
            global simulateData
            global fitData
            global saveResults
            global abortNow
            
            Carrier_Starting
            % Get the user specified MAT File
            DataToAnalyze1 = app.TRPLMATfilePath.Value;
            DataToAnalyze2 = app.PLMATfilePath.Value;
            
           
            noPLData = app.NoPLDataCheckBox.Value;
            % MAT File should contain the following
            %   Time_Axis:         Time axis of acquisition, column vector in units of us
            %   Temperatures:      Temperatures for each acquisition, row vector in units of K
            %   Excitations:       Excitations for each acquisition, col vector, units of deg (TRPL) or V (TMR)
            %   AdaptiveParameter: Adaptable parameter that can be swept through
            %   Adaptives:         Values of the adaptable parameter for each acquisition, col vector   
            %   RawData:           Vector containing raw data, size = (Time_Axis,Temperatures,Excitations,Adaptives)
            %   SampleName:        Sample name and identifying information
            %   TRPL_or_TMR:       String specifying if loaded data is TRPL or TMR acquisition
            %   Decay:             Values of the decay data with the baseline subtracted
            %   Decay_WithinFitWindow: Fit window chosen for the decay data
            %   t0_vTemp:          lifetime calculated to extract the data 
            
            
          
            try
                load( DataToAnalyze1 , 'Time_Axis' , 'Temperatures' , 'Excitations' ,...
                    'SampleName' , 'TRPL_or_TMR' , 'AdaptiveParameter' ,...
                    'Adaptives' , 'AdaptiveUnits',  't0_vTemp');
            catch
                ErrorMessage = 'Problem trying to load TRPL MAT file at specified path';
                app.ErrorReport.Value = ErrorMessage;
                error(ErrorMessage)
            end
            
          
            
            
            
          
          

            
            
            % Create a color wheel to overlay temperature-dependent data, different colors with each excitation level
            ColorWheel = [ 'bo-' ; 'go-' ; 'ro-' ; 'ko-' ; 'co-' ; 'yo-' ];
            
            % Create a color wheel to overlay temperature-dependent data, different styles with each
            % recombination component
            ColorWheel2 = [ "k--" ; "k-." ; "k:" ; "k-"; "r-."];
  
            
            adaptive_Ind = 1;
            
            %Plot as a function of excitations
            for exciteInd = 1:length(Excitations)
                plot(app.UIAxes,Temperatures, t0_vTemp(adaptive_Ind,:,exciteInd), ColorWheel(exciteInd,:),...
                    'Displayname', string(Excitations(exciteInd)))
                hold(app.UIAxes,'on')
            end
            
            hold(app.UIAxes,'off')
            
            
            % Initialize parameters that will be used in the simulation
            % The database from the CSV will be read and temperature probes
            % are initialized here too
            materialsDataBase = readtable('BandParameters.csv');
            
            % have a temperature simulation from 50-330 K in steps of 1 K
            beginTemp = 50;
            endTemp = 330;
            stepNumber = endTemp - beginTemp +1;
            Tprobe = linspace(beginTemp,endTemp,stepNumber)     ;
            tauAug = [];
            tauRad = [];
            tauSRH = [];
            totalLifetime = [];
         
            
            
            
            
            simulateData = 0;
            fitData = 0;
            saveResults = 0;
            abortNow = 0;
           
           
            % Reset WaitingForUser to 1 to hold in while-loop
            waitingForUser = 1;
                               
            % WhileBreaker variable increments to time out the while-loop
            whileBreaker = 0;
                                
            % Define while-loop time increment
            whileIterationTime = 0.05; % 50 ms iterations
                                
            % Define while-loop total time limit to break loop
            whileTimeLimit = 360; % 360 seconds = 5 minutes
            whileIterationLimit = whileTimeLimit/whileIterationTime;
            
            % Operations started, set status to active
            Carrier_Activate;
            while(waitingForUser)
                
              % Wait the while-loop time incre
                pause(whileIterationTime) 
                                    
              % Increment WhileBreaker
                whileBreaker = whileBreaker + 1;
                
                
                
                
              % Simulate the data from the parameters box
                if(simulateData)
                    material = '';
                    bandParams = [];
                    materialDopingType = '';
                    dopingDensity = nan;
                    defectDensity = nan;
                    defectLevel  = nan;
                    DataToAnalyze2 = '';
                    Nc  = nan;
                    Nv  = nan;
                    ni  = nan;
                    G   = nan;
                    phi = nan;
                    
                    
                    
                    
                 %   materialsDataBase = readtable('BandParameters.csv');  
                    material = app.MaterialDropDown.Value;
                    bandParams = nan*ones(9,1);
                   
           
                    switch(material) 
                     case 'InAs'
                        bandParams =materialsDataBase.InAs;
       
                     case 'GaAs'
                        bandParams = materialsDataBase.GaAs;
                     case 'InSb'
                        bandParams = materialsDataBase.InSb;
                     case 'InAsSb'
                         % use InAs band parameters
                        bandParams = materialsDataBase.InAs;
                     otherwise
                        bandParams = materials.DataBase.InAs;
                
                    end
                    % Band Parameters 
                    %row1     =  Band gap in eV (0 K) perhaps
                    %row2     =  electron effective mass in units of free electron mass
                    %row3     =  gamma1 inverse effective mass parameter 1
                    %row4     =  gamma2 inverse effective mass parameter 2
                    %row5     =  gamma3 inverse effective mass parameter 3
                    %row6     =  alpha 0 K band gap using Varshni equation
                    %row7     =  beta  slope of the Varshni equation
                    %row8     =  einf static dielectric constant
                    %row9     =  F1F2  Auger overlap integral parameter
            
                    type = 'i';
                    materialDopingType = app.DopingDropDown.Value;
                      % Check the doping Type
                    switch(materialDopingType)
                     case('Intrinsic')
                        type = 'i';
                     case('p-type')
                        type ='p';
                     case('n-type')
                        type = 'n';
                    
                    end              
                    %Doping density input will be in units of 1e15 cm^-3
                    dopingDensity = app.DopingDensityEditField.Value;
                    
                    %Defect density in units of m^-1
                    defectDensity = app.DefectDensity1EditField.Value;
            
                    % Units of meV, so change to eV The negative shows that
                    % the conduction band edge is at zero eV
                    defectLevel = -app.DefectLevel1EditField.Value/1000;
                    
                    % Bloch overlap integral value
                    blochOverlap = app.BlochoverlapEditField.Value;
                    
                    % For the case of intentional doping
                    defectDensity2 = app.DefectDensity2EditField.Value;
                    defectLevel2 = -app.DefectLevel2EditField.Value/1000;
                    
                    
                    
                    % extract the band parameters and use the Varshni equation
                    % with given parameters from the Vurgaftman review article
                    % See the CSV named bandParameters.csv to look at the
                    % input parameters
                    [meStar, mhStar, eg, nRefractive,...
                kExtinction, photonEnergyum, varshniEinstein, einf, f1f2] = extractParameters(bandParams, material);
            
                    
                
                
                
                
                    % If there is PL data, we will use the Bose Einstein single
                    % oscillator equation 
                    % The exponential modified gaussian fit method will be used, which is
                    % close to the first derivative method of extracting the band
                    % gap
                    
                    
                    % Load PL fitted data and use the bose einstein parameters for
                    % bandgap as a function of temperature
                    noPLData = app.NoPLDataCheckBox.Value;
                    DataToAnalyze2 = app.PLMATfilePath.Value;
                    if(~noPLData)
                        try
                            load(DataToAnalyze2, 'boseParameters')
                             varshniEinstein = @(T) boseEinstein(T, boseParameters(1,1), boseParameters(2,1), ...
                             boseParameters(3,1));
                             eg = boseParameters(1,1);
                   
                        catch 
                            ErrorMessage = 'Problem trying to load PL MAT file: does file exist and contain "boseParameters?"';
                            app.ErrorReport.Value = ErrorMessage;
                            error(ErrorMessage);
                    
                        end
                    end
                  
                    % Band gap at 300 K in eV
                    eg300 = varshniEinstein(300) / 1000;
            
                    % Temperature dependent band gap in units of eV
                    egTemp = varshniEinstein(Tprobe) ./ 1000;
           
                    % Conduction band edge energy will be the reference
                    % point
                    % 0K
                    conductionEdge = 0;
            
                    % Set the valence band edge having the bose einstein temperature
                    % dependence IT has a negative value
                    valenceEdge = conductionEdge - egTemp;
                
                    %get the active region thickness and convert it to nm
                    t = app.ActiveRegionThicknessEditField.Value * 1000;
                    n = app.CaplayerrefractiveindexatbandgapenergyEditField.Value;
                    
                    [Nc, Nv, ni, G, phi] = ...
                        calculateKeyParameters(Tprobe, meStar, mhStar, photonEnergyum,...
                        nRefractive,kExtinction, eg300, egTemp, 1200, n, t,einf );
                    
                   intentionalDoping = app.IntentionaldopingCheckBox.Value; 
                        
                   
                   if(intentionalDoping)
                         % calculate the fitted lifetimes
                    [totalLifetime, tauRad, tauSRH,tauSRH2, tauAug] = calculateLifetimes2(Tprobe, type,meStar,...
                    mhStar, egTemp, valenceEdge,conductionEdge, einf,...
                    Nc, Nv, ni, G, phi, blochOverlap, dopingDensity, defectLevel, defectDensity, defectLevel2, defectDensity2);
            
                       
                   else
                       
                       % calculate the fitted lifetimes
                    [totalLifetime, tauRad, tauSRH, tauAug] = calculateLifetimes(Tprobe, type,meStar,...
                    mhStar, egTemp, valenceEdge,conductionEdge, einf,...
                    Nc, Nv, ni, G, phi, blochOverlap, dopingDensity, defectLevel, defectDensity);
            
                       
                       
                   end
                   
                    
                    %should we plot in a new Figure?
                    plotresults = app.ShowResultsCheckBox.Value;

                    
                    if(plotresults)
                        %PLot the lifetime components in a new plot
                        
                        figure
                        hold on
                        for exciteInd = 1:length(Excitations)
                            plot(Temperatures, t0_vTemp(adaptive_Ind,:,exciteInd), ColorWheel(exciteInd,:),...
                                'Displayname', string(Excitations(exciteInd)))
                        end
                        plot(Tprobe,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                        plot(Tprobe,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                        
                        if(intentionalDoping)
                        plot(Tprobe,tauSRH2,ColorWheel2(5,:), 'Displayname', 'SRH2')
                        end
                        
                        plot(Tprobe,tauAug,ColorWheel2(3,:), 'Displayname', 'Auger')
                        plot(Tprobe,totalLifetime, ColorWheel2(4,:),'Displayname','Total')
                        legend()
                        set(gca,'YScale', 'log')
                        hold off
                        ylim([1e-1, 1e3])
                        
                    else
                
                        % Plot as a function of excitations
                        for exciteInd = 1:length(Excitations)
                            plot(app.UIAxes,Temperatures, t0_vTemp(adaptive_Ind,:,exciteInd), ColorWheel(exciteInd,:),...
                                'Displayname', string(Excitations(exciteInd)))
                            hold(app.UIAxes,'on')
                        end
                    
                        %Plot the lifetime components in the app window
                        plot(app.UIAxes,Tprobe,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                        plot(app.UIAxes,Tprobe,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                        if(intentionalDoping)
                        plot(app.UIAxes,Tprobe,tauSRH2,ColorWheel2(5,:), 'Displayname', 'SRH2')
                        end 
                        plot(app.UIAxes,Tprobe,tauAug,ColorWheel2(3,:), 'Displayname', 'Auger')
                        plot(app.UIAxes,Tprobe,totalLifetime, ColorWheel2(4,:),'Displayname','Total')
                        legend(app.UIAxes)
                        set(app.UIAxes, 'YScale', 'log')
                        hold(app.UIAxes, 'off')
                        ylim(app.UIAxes,[1e-1, 1e3])
                        
                    end
                    
                    
                    % Go back to waiting mode
                    simulateData = 0;
                    
                    %Restart while Breaker
                    whileBreaker = 0;

                    
                    %It will essentially do the same thing but first fit
                    %and then plot
                elseif(fitData)
                    % re initialize the data, because user may go straight
                    % to fitting routine
                    material = '';
                    bandParams = [];
                    materialDopingType = '';
                    dopingDensity = nan;
                    defectDensity = nan;
                    defectLevel  = nan;
                    DataToAnalyze2 = '';
                    Nc  = nan;
                    Nv  = nan;
                    ni  = nan;
                    G   = nan;
                    phi = nan;
                  %  materialsDataBase = readtable('BandParameters.csv');  
                    material = app.MaterialDropDown.Value;
                   
                    %row1     =  Band gap in eV (0 K) perhaps
                    %row2     =  electron effective mass in units of free electron mass
                    %row3     =  gamma1 inverse effective mass parameter 1
                    %row4     =  gamma2 inverse effective mass parameter 2
                    %row5     =  gamma3 inverse effective mass parameter 3
           
                    switch(material) 
                     case 'InAs'
                        bandParams =materialsDataBase.InAs;
       
                     case 'GaAs'
                        bandParams = materialsDataBase.GaAs;
                     case 'InSb'
                        bandParams = materialsDataBase.InSb;
                     case 'InAsSb'
                         % use InAs band parameters
                        bandParams = materialsDataBase.InAs;
                     otherwise
                        bandParams = materials.DataBase.InAs;
                
                    end
            
                    type = 'i';
                    materialDopingType = app.DopingDropDown.Value;
                      % Check the doping Type
                    switch(materialDopingType)
                     case('Intrinsic')
                        type = 'i';
                     case('p-type')
                        type ='p';
                     case('n-type')
                        type = 'n';
                    
                    end    
                    
                    
                    
                    % in units of 1e15 per cm^3 we will convert to per cm^3
                    % in the functions (fitting a very large number seems
                    % to cause it to be insensitive)
                    dopingDensity = app.DopingDensityEditField.Value;
            
                    % in units of per m
                    defectDensity = app.DefectDensity1EditField.Value;
            
                    % Units of meV, so change to eV The negative shows that
                    % the conduction band edge is the zero energy level
                    defectLevel = -app.DefectLevel1EditField.Value/1000;
                    
                    
                    % For the case of intentional doping
                    defectDensity2 = app.DefectDensity2EditField.Value;
                    defectLevel2 = -app.DefectLevel2EditField.Value/1000;
            
            
                    intentionalDoping = app.IntentionaldopingCheckBox.Value;
                  
                    % extract the band parameters and use the Varshni equation
                    % with given parameters from the Vurgaftman review article
                    [meStar, mhStar, eg, nRefractive,...
                kExtinction, photonEnergyum, varshniEinstein, einf, f1f2] = extractParameters(bandParams, material);
            
                
                    % If there is PL data, we will use the Bose Einstein single
                    % oscillator equation 
                    % The exponential modified gaussian fit method will be used, which is
                    % close in agreement to the first derivative method of extracting the band
                    % gap----SEE the PL_Analyzer app in PL_AnalyzerCurveFittingToolboxV2
                    noPLData = app.NoPLDataCheckBox.Value;
                    if(~noPLData)
                        varshniEinstein = @(T) boseEinstein(T, boseParameters(1,1), boseParameters(2,1), ...
                        boseParameters(3,1));
                        eg = boseParameters(1,1);
                    end
             
                    
            
                    % Band gap at 300 K in eV
                    eg300 = varshniEinstein(300) / 1000;
            
                    % Temperature dependent band gap in units of eV
                    egTemp = varshniEinstein(Tprobe) ./ 1000;
                    egTemp2 = varshniEinstein(Temperatures)./ 1000;
           
                    % Conduction band edge energy will be equal to the band gap at
                    % 0K
                    conductionEdge = 0;
            
                    % Set the valence band edge having the bose einstein temperature
                    % dependence
                    valenceEdge =  conductionEdge - egTemp;
                    valenceEdge2 = conductionEdge - egTemp2;
                
                    %get the active region thickness and convert it to nm
                    t = app.ActiveRegionThicknessEditField.Value * 1000;
                    n = app.CaplayerrefractiveindexatbandgapenergyEditField.Value;
                    
                    %The first excitation will be used to fit the lifetime
                    % 
                    y = t0_vTemp(adaptive_Ind,:,1);
                    
                     [Nc2, Nv2, ni2, G2, phi] = ...
                        calculateKeyParameters(Temperatures, meStar, mhStar, photonEnergyum,...
                        nRefractive,kExtinction, eg300, egTemp2, 1200, n, t, einf);
                    
                    % Take the user input auger overlap value which is 0.15
                    blochOverlap = app.BlochoverlapEditField.Value;
             
                    augeroverlapfit = app.FitCheckBox.Value;
                    fix = app.FixCheckBox.Value;
                    if(intentionalDoping)
                    [fobject, gof, output] = fitLifetimes(Temperatures, y,augeroverlapfit,fix,type,...
                    meStar,mhStar, egTemp2, valenceEdge2, conductionEdge, einf,...
                    Nc2, Nv2, ni2, G2, phi, blochOverlap, dopingDensity, defectLevel, defectDensity,defectLevel2, defectDensity2);
                        defectLevel2 = fobject.s;
                        app.DefectLevel2EditField.Value = -defectLevel2*1000;
                        
                        defectDensity2 = fobject.t;
                        app.DefectDensity2EditField.Value = defectDensity2;
                    else
                    
                    [fobject, gof, output] = fitLifetimes(Temperatures, y,augeroverlapfit,fix,type,...
                    meStar,mhStar, egTemp2, valenceEdge2, conductionEdge, einf,...
                    Nc2, Nv2, ni2, G2, phi, blochOverlap, dopingDensity, defectLevel, defectDensity);
                    
                    end
                    
                    
                    dopingDensity = fobject.p;
                    app.DopingDensityEditField.Value = dopingDensity;
                    
                    defectLevel = fobject.q;
                    app.DefectLevel1EditField.Value = -defectLevel * 1000;
                    
                    defectDensity = fobject.r;
                    app.DefectDensity1EditField.Value = defectDensity;
                    
                    blochOverlap = fobject.o;
                    app.BlochoverlapEditField.Value = blochOverlap;
                    
                    
                    
                    
                    
                    
                    
                     [Nc, Nv, ni, G, phi] = ...
                        calculateKeyParameters(Tprobe, meStar, mhStar, photonEnergyum,...
                        nRefractive,kExtinction, eg300, egTemp, 1200, n, t, einf);
                    
                        
                        
                    if (intentionalDoping)
                        [totalLifetime, tauRad, tauSRH, tauSRH2, tauAug] = calculateLifetimes2(Tprobe, type,meStar,...
                    mhStar, egTemp, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlap, dopingDensity, defectLevel, defectDensity, defectLevel2, defectDensity2);
            
                    else
                        % calculate the lifetimes
                    [totalLifetime, tauRad, tauSRH, tauAug] = calculateLifetimes(Tprobe, type,meStar,...
                    mhStar, egTemp, valenceEdge, conductionEdge, einf,Nc, Nv, ni, G, phi,...
                    blochOverlap, dopingDensity, defectLevel, defectDensity);
            
                        
                    end
                    
                
                    % Plot as a function of excitations
                    for exciteInd = 1:length(Excitations)
                        plot(app.UIAxes,Temperatures, t0_vTemp(adaptive_Ind,:,exciteInd), ColorWheel(exciteInd,:),...
                             'Displayname', string(Excitations(exciteInd)))
                        hold(app.UIAxes,'on')
                    end
                    plot(app.UIAxes,Tprobe,tauRad,ColorWheel2(1,:),'Displayname', 'Radiative')
                    
                    plot(app.UIAxes,Tprobe,tauSRH,ColorWheel2(2,:), 'Displayname', 'SRH')
                    if(intentionalDoping)
                    plot(app.UIAxes,Tprobe,tauSRH2,ColorWheel2(5,:), 'Displayname', 'SRH2')

                    end
                    
                    plot(app.UIAxes,Tprobe,tauAug,ColorWheel2(3,:), 'Displayname', 'Auger')
                    plot(app.UIAxes,Tprobe,totalLifetime, ColorWheel2(4,:),'Displayname','Total')
                    legend(app.UIAxes)
                    set(app.UIAxes, 'YScale', 'log')
                    
                    hold(app.UIAxes, 'off')
                    ylim(app.UIAxes,[1e-1, 1e3])
                    
                    
                    
                   
                    
                    
                    
                    %Reset the fitting routine
                    fitData = 0;
                    
                    %Reset the wait time to exit/abort
                    whileBreaker = 0;
                    
                    
                elseif(saveResults)
                    saveTable = table();
                    intentionalDoping =app.IntentionaldopingCheckBox.Value;
                    try
                        saveTable = table(Tprobe',tauRad,tauSRH,tauAug,totalLifetime);
                    catch 
                        errorMessage = 'Problem Creating a table to save lifetime components. Try to simulate data';
                        app.ErrorReport.Value = errorMessage;
                        saveResults = 0;
                        continue;
                        
                    end
                    LT_Writing;
                    if(intentionalDoping)
                         saveTable = table(Tprobe',tauRad,tauSRH,tauSRH2,tauAug,totalLifetime);
                         saveTable.Properties.VariableNames = {'Temp_K', 'tauRad_us', 'tauSRH_us','tauSRH2_us' 'tauAug_us','Total_us'};
                         filename = string(SampleName) + "_LifetimeComponents.txt";
                         writetable(saveTable, filename);
                    
                   
                         saveTable2 = table(dopingDensity,app.DefectLevel1EditField.Value, defectDensity,...
                             app.DefectLevel2EditField.Value, app.DefectDensity2EditField.value);
                         
                         saveTable2.Properties.VariableNames = {'DopingDensity_1e15_cm3', 'DefectLevel_meV', 'DefectDensity_m1',...
                             'DefectLevel2_meV', 'DefectDensity2_m1'};
                         filename2 = string(SampleName) + '_InputParameters.txt';
                         writetable(saveTable2,filename2);
                    
                    else
                         saveTable.Properties.VariableNames = {'Temp_K', 'tauRad_us', 'tauSRH_us', 'tauAug_us','Total_us'};
                         filename = string(SampleName) + "_LifetimeComponents.txt";
                         writetable(saveTable, filename);
                    
                   
                         saveTable2 = table(dopingDensity,app.DefectLevel1EditField.Value, defectDensity);
                         saveTable2.Properties.VariableNames = {'DopingDensity_1e15_cm3', 'DefectLevel_meV', 'defectDensity_m1'};
                         filename2 = string(SampleName) + '_InputParameters.txt';
                         writetable(saveTable2,filename2);
                    
                    
                    end
                    Carrier_Activate;
                    %Reset the save routine and leave app in wait mode
                    saveResults = 0;
                    
                    
                    
                    whileBreaker = 0;
                    
            
                elseif(abortNow || (whileBreaker > whileIterationLimit))
                    % Operations aborted, reset tool
                      LT_HardReset
                                      
                    % Reset clears the decay plot, also clear the temperature dependence plot
                      scatter( app.UIAxes , 1 , 1 , 'x' ) 
                      abortNow = 0;
                    % Abort if user presses button or if whileIterationLimit time-out reached                                
                      return
                    
                    
                    
                    
                    
                    
                    
                end
               
             end
                

                            
        
           
           
            
        end

        % Button pushed function: AbortButton
        function AbortCalculation(app, event)
            
            % Define global variable AbortNow = 1 to enable analysis to repeat
            global abortNow
            abortNow = 1;
                  
            
        end

        % Button pushed function: SimulateButton
        function simulateLifetimes(app, event)
              
          global simulateData
          simulateData = 1;
        
        end

        % Button pushed function: FitButton
        function fitLifetimes(app, event)
            global fitData
            fitData = 1;
        end

        % Button pushed function: SaveButton
        function saveModel(app, event)
            global saveResults
            saveResults = 1;
        
        end

        % Value changed function: IntentionaldopingCheckBox
        function Doping(app, event)
            value = app.IntentionaldopingCheckBox.Value;
            if(value)
                app.DefectDensity2EditField.Editable = 'on';
                app.DefectDensity2EditField.Enable = 'on';
                
                app.DefectLevel2EditField.Editable = 'on';
                app.DefectLevel2EditField.Enable ='on';
                
                app.FixCheckBox.Value = true;
                app.DefectDensity1EditField.Editable = 'off';    
                app.DefectLevel1EditField.Editable = 'off';
                app.FixCheckBox.Enable ='on';
                
            else
                % Make the editfield seem like it won't be used
                app.DefectDensity2EditField.Editable = 'off';
                app.DefectDensity2EditField.Enable = 'off';
      
                %Make the editfield seem like it won't be used
                app.DefectLevel2EditField.Editable = 'off';
                app.DefectLevel2EditField.Enable = 'off';
                
                % The first defect Level will not be used to
                app.FixCheckBox.Value = false;
                app.DefectDensity1EditField.Editable = 'on';
                app.DefectLevel1EditField.Editable = 'on';
                app.FixCheckBox.Enable = 'off';
                
                
                
            end
        end

        % Value changed function: FixCheckBox
        function FixDefect1(app, event)
            value = app.FixCheckBox.Value;
            if(value)
                app.DefectDensity1EditField.Editable = 'off';    
                app.DefectLevel1EditField.Editable = 'off';
            else 
                app.DefectDensity1EditField.Editable = 'on';
                app.DefectLevel1EditField.Editable = 'on'; 
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 698 559];
            app.UIFigure.Name = 'UI Figure';

            % Create DataFittingStatusEditFieldLabel
            app.DataFittingStatusEditFieldLabel = uilabel(app.UIFigure);
            app.DataFittingStatusEditFieldLabel.HorizontalAlignment = 'right';
            app.DataFittingStatusEditFieldLabel.FontSize = 16;
            app.DataFittingStatusEditFieldLabel.FontWeight = 'bold';
            app.DataFittingStatusEditFieldLabel.Position = [70 524 152 22];
            app.DataFittingStatusEditFieldLabel.Text = 'Data Fitting Status:';

            % Create ToolStatus
            app.ToolStatus = uieditfield(app.UIFigure, 'text');
            app.ToolStatus.Editable = 'off';
            app.ToolStatus.HorizontalAlignment = 'center';
            app.ToolStatus.FontSize = 14;
            app.ToolStatus.FontWeight = 'bold';
            app.ToolStatus.FontColor = [1 0 0];
            app.ToolStatus.Position = [237 524 134 22];
            app.ToolStatus.Value = 'INACTIVE';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Temperature Dependence')
            xlabel(app.UIAxes, 'Temperature (K)')
            ylabel(app.UIAxes, 'Characteristic Slope (us)')
            app.UIAxes.PlotBoxAspectRatio = [1.09122807017544 1 1];
            app.UIAxes.Box = 'on';
            app.UIAxes.Position = [9 188 362 328];

            % Create InputparametersandusercontrolPanel
            app.InputparametersandusercontrolPanel = uipanel(app.UIFigure);
            app.InputparametersandusercontrolPanel.Title = 'Input parameters and user control';
            app.InputparametersandusercontrolPanel.Position = [376 182 304 343];

            % Create MaterialDropDownLabel
            app.MaterialDropDownLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.MaterialDropDownLabel.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel.Position = [11 297 48 22];
            app.MaterialDropDownLabel.Text = 'Material';

            % Create MaterialDropDown
            app.MaterialDropDown = uidropdown(app.InputparametersandusercontrolPanel);
            app.MaterialDropDown.Items = {'InAs', 'GaAs', 'InAsSb', 'InSb'};
            app.MaterialDropDown.Position = [74 297 100 22];
            app.MaterialDropDown.Value = 'InAs';

            % Create IntentionaldopingCheckBox
            app.IntentionaldopingCheckBox = uicheckbox(app.InputparametersandusercontrolPanel);
            app.IntentionaldopingCheckBox.ValueChangedFcn = createCallbackFcn(app, @Doping, true);
            app.IntentionaldopingCheckBox.Text = 'Intentional doping?';
            app.IntentionaldopingCheckBox.Position = [179 297 123 22];

            % Create DopingTypeDropDownLabel
            app.DopingTypeDropDownLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DopingTypeDropDownLabel.HorizontalAlignment = 'right';
            app.DopingTypeDropDownLabel.Position = [3 266 73 22];
            app.DopingTypeDropDownLabel.Text = 'Doping Type';

            % Create DopingDropDown
            app.DopingDropDown = uidropdown(app.InputparametersandusercontrolPanel);
            app.DopingDropDown.Items = {'Intrinsic', 'n-type', 'p-type', 'Option 4'};
            app.DopingDropDown.Position = [91 266 138 22];
            app.DopingDropDown.Value = 'n-type';

            % Create DopingDensityEditFieldLabel
            app.DopingDensityEditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DopingDensityEditFieldLabel.HorizontalAlignment = 'right';
            app.DopingDensityEditFieldLabel.Position = [4 235 86 22];
            app.DopingDensityEditFieldLabel.Text = 'Doping Density';

            % Create DopingDensityEditField
            app.DopingDensityEditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.DopingDensityEditField.Position = [105 235 81 22];
            app.DopingDensityEditField.Value = 1.7;

            % Create X1e15cm3Label
            app.X1e15cm3Label = uilabel(app.InputparametersandusercontrolPanel);
            app.X1e15cm3Label.Position = [199 235 78 22];
            app.X1e15cm3Label.Text = 'X1e15 cm(-3)';

            % Create DefectLevel1EditFieldLabel
            app.DefectLevel1EditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DefectLevel1EditFieldLabel.HorizontalAlignment = 'right';
            app.DefectLevel1EditFieldLabel.Position = [4 206 82 22];
            app.DefectLevel1EditFieldLabel.Text = 'Defect Level 1';

            % Create DefectLevel1EditField
            app.DefectLevel1EditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.DefectLevel1EditField.Position = [116 206 81 22];
            app.DefectLevel1EditField.Value = 67;

            % Create meVLabel_7
            app.meVLabel_7 = uilabel(app.InputparametersandusercontrolPanel);
            app.meVLabel_7.Position = [210 206 30 22];
            app.meVLabel_7.Text = 'meV';

            % Create FixCheckBox
            app.FixCheckBox = uicheckbox(app.InputparametersandusercontrolPanel);
            app.FixCheckBox.ValueChangedFcn = createCallbackFcn(app, @FixDefect1, true);
            app.FixCheckBox.Enable = 'off';
            app.FixCheckBox.Text = 'Fix?';
            app.FixCheckBox.Position = [243 206 45 22];

            % Create DefectDensity1EditFieldLabel
            app.DefectDensity1EditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DefectDensity1EditFieldLabel.HorizontalAlignment = 'right';
            app.DefectDensity1EditFieldLabel.Position = [4 174 94 22];
            app.DefectDensity1EditFieldLabel.Text = 'Defect Density 1';

            % Create DefectDensity1EditField
            app.DefectDensity1EditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.DefectDensity1EditField.Position = [116 174 81 22];
            app.DefectDensity1EditField.Value = 0.7;

            % Create m1Label_7
            app.m1Label_7 = uilabel(app.InputparametersandusercontrolPanel);
            app.m1Label_7.Position = [210 174 34 22];
            app.m1Label_7.Text = 'm(-1)';

            % Create DefectLevel2EditFieldLabel
            app.DefectLevel2EditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DefectLevel2EditFieldLabel.HorizontalAlignment = 'right';
            app.DefectLevel2EditFieldLabel.Position = [10 143 82 22];
            app.DefectLevel2EditFieldLabel.Text = 'Defect Level 2';

            % Create DefectLevel2EditField
            app.DefectLevel2EditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.DefectLevel2EditField.Editable = 'off';
            app.DefectLevel2EditField.Enable = 'off';
            app.DefectLevel2EditField.Position = [115 143 87 22];
            app.DefectLevel2EditField.Value = 70;

            % Create meVLabel_8
            app.meVLabel_8 = uilabel(app.InputparametersandusercontrolPanel);
            app.meVLabel_8.Position = [214 143 30 22];
            app.meVLabel_8.Text = 'meV';

            % Create DefectDensity2EditFieldLabel
            app.DefectDensity2EditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.DefectDensity2EditFieldLabel.HorizontalAlignment = 'right';
            app.DefectDensity2EditFieldLabel.Position = [6 113 94 22];
            app.DefectDensity2EditFieldLabel.Text = 'Defect Density 2';

            % Create DefectDensity2EditField
            app.DefectDensity2EditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.DefectDensity2EditField.Editable = 'off';
            app.DefectDensity2EditField.Enable = 'off';
            app.DefectDensity2EditField.Position = [115 113 100 22];
            app.DefectDensity2EditField.Value = 50;

            % Create m1Label_8
            app.m1Label_8 = uilabel(app.InputparametersandusercontrolPanel);
            app.m1Label_8.Position = [221 113 34 22];
            app.m1Label_8.Text = 'm(-1)';

            % Create BlochoverlapEditFieldLabel
            app.BlochoverlapEditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.BlochoverlapEditFieldLabel.HorizontalAlignment = 'right';
            app.BlochoverlapEditFieldLabel.Position = [14 83 78 22];
            app.BlochoverlapEditFieldLabel.Text = 'Bloch overlap';

            % Create BlochoverlapEditField
            app.BlochoverlapEditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.BlochoverlapEditField.Position = [112 83 82 22];
            app.BlochoverlapEditField.Value = 0.15;

            % Create FitCheckBox
            app.FitCheckBox = uicheckbox(app.InputparametersandusercontrolPanel);
            app.FitCheckBox.Text = 'Fit?';
            app.FitCheckBox.Position = [219 83 42 22];

            % Create ActiveRegionThicknessEditFieldLabel
            app.ActiveRegionThicknessEditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.ActiveRegionThicknessEditFieldLabel.HorizontalAlignment = 'right';
            app.ActiveRegionThicknessEditFieldLabel.Position = [4 47 137 22];
            app.ActiveRegionThicknessEditFieldLabel.Text = 'Active Region Thickness';

            % Create ActiveRegionThicknessEditField
            app.ActiveRegionThicknessEditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.ActiveRegionThicknessEditField.Position = [156 47 61 22];
            app.ActiveRegionThicknessEditField.Value = 2;

            % Create umLabel
            app.umLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.umLabel.Position = [229 47 30 22];
            app.umLabel.Text = '(um)';

            % Create CaplayerrefractiveindexatbandgapenergyEditFieldLabel
            app.CaplayerrefractiveindexatbandgapenergyEditFieldLabel = uilabel(app.InputparametersandusercontrolPanel);
            app.CaplayerrefractiveindexatbandgapenergyEditFieldLabel.HorizontalAlignment = 'right';
            app.CaplayerrefractiveindexatbandgapenergyEditFieldLabel.Position = [2 17 247 22];
            app.CaplayerrefractiveindexatbandgapenergyEditFieldLabel.Text = 'Cap layer refractive index at bandgap energy';

            % Create CaplayerrefractiveindexatbandgapenergyEditField
            app.CaplayerrefractiveindexatbandgapenergyEditField = uieditfield(app.InputparametersandusercontrolPanel, 'numeric');
            app.CaplayerrefractiveindexatbandgapenergyEditField.Position = [261 17 38 22];
            app.CaplayerrefractiveindexatbandgapenergyEditField.Value = 3.7;

            % Create DataFittingPanel
            app.DataFittingPanel = uipanel(app.UIFigure);
            app.DataFittingPanel.Title = 'Data Fitting';
            app.DataFittingPanel.FontWeight = 'bold';
            app.DataFittingPanel.Position = [9 46 671 128];

            % Create TRPLFittedresultsMATfileLabel
            app.TRPLFittedresultsMATfileLabel = uilabel(app.DataFittingPanel);
            app.TRPLFittedresultsMATfileLabel.HorizontalAlignment = 'right';
            app.TRPLFittedresultsMATfileLabel.Position = [10 78 158 22];
            app.TRPLFittedresultsMATfileLabel.Text = 'TRPL Fitted results MAT file:';

            % Create TRPLMATfilePath
            app.TRPLMATfilePath = uieditfield(app.DataFittingPanel, 'text');
            app.TRPLMATfilePath.Position = [173 78 157 22];
            app.TRPLMATfilePath.Value = 'L19-091_FittedResults.mat';

            % Create PLFittedresultsMATfileEditFieldLabel
            app.PLFittedresultsMATfileEditFieldLabel = uilabel(app.DataFittingPanel);
            app.PLFittedresultsMATfileEditFieldLabel.HorizontalAlignment = 'right';
            app.PLFittedresultsMATfileEditFieldLabel.Position = [337 78 142 22];
            app.PLFittedresultsMATfileEditFieldLabel.Text = 'PL Fitted results MAT file:';

            % Create PLMATfilePath
            app.PLMATfilePath = uieditfield(app.DataFittingPanel, 'text');
            app.PLMATfilePath.Position = [490 78 157 22];
            app.PLMATfilePath.Value = 'L19-091_FittedResultsPL.mat';

            % Create StartButton
            app.StartButton = uibutton(app.DataFittingPanel, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @Initialize, true);
            app.StartButton.FontWeight = 'bold';
            app.StartButton.FontColor = [0.3922 0.8314 0.0745];
            app.StartButton.Position = [20 42 71 22];
            app.StartButton.Text = 'Start';

            % Create SimulateButton
            app.SimulateButton = uibutton(app.DataFittingPanel, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @simulateLifetimes, true);
            app.SimulateButton.Position = [100 42 74 22];
            app.SimulateButton.Text = 'Simulate';

            % Create FitButton
            app.FitButton = uibutton(app.DataFittingPanel, 'push');
            app.FitButton.ButtonPushedFcn = createCallbackFcn(app, @fitLifetimes, true);
            app.FitButton.Position = [183 42 74 22];
            app.FitButton.Text = 'Fit';

            % Create SaveButton
            app.SaveButton = uibutton(app.DataFittingPanel, 'push');
            app.SaveButton.ButtonPushedFcn = createCallbackFcn(app, @saveModel, true);
            app.SaveButton.Position = [264 42 74 22];
            app.SaveButton.Text = 'Save';

            % Create AbortButton
            app.AbortButton = uibutton(app.DataFittingPanel, 'push');
            app.AbortButton.ButtonPushedFcn = createCallbackFcn(app, @AbortCalculation, true);
            app.AbortButton.FontColor = [1 0 0];
            app.AbortButton.Position = [350 42 85 22];
            app.AbortButton.Text = 'Abort';

            % Create NoPLDataCheckBox
            app.NoPLDataCheckBox = uicheckbox(app.DataFittingPanel);
            app.NoPLDataCheckBox.Text = 'No PL Data?';
            app.NoPLDataCheckBox.Position = [525 53 91 22];

            % Create OutputfilesuffixLabel
            app.OutputfilesuffixLabel = uilabel(app.DataFittingPanel);
            app.OutputfilesuffixLabel.HorizontalAlignment = 'right';
            app.OutputfilesuffixLabel.Position = [10 9 95 22];
            app.OutputfilesuffixLabel.Text = 'Output file suffix:';

            % Create UserOutputFileSuffix
            app.UserOutputFileSuffix = uieditfield(app.DataFittingPanel, 'text');
            app.UserOutputFileSuffix.Position = [108 9 118 22];

            % Create ShowResultsCheckBox
            app.ShowResultsCheckBox = uicheckbox(app.DataFittingPanel);
            app.ShowResultsCheckBox.Text = 'Show Results?';
            app.ShowResultsCheckBox.Position = [244 9 102 22];

            % Create ErrorEditFieldLabel
            app.ErrorEditFieldLabel = uilabel(app.UIFigure);
            app.ErrorEditFieldLabel.HorizontalAlignment = 'right';
            app.ErrorEditFieldLabel.FontWeight = 'bold';
            app.ErrorEditFieldLabel.Position = [9 15 39 22];
            app.ErrorEditFieldLabel.Text = 'Error:';

            % Create ErrorReport
            app.ErrorReport = uieditfield(app.UIFigure, 'text');
            app.ErrorReport.Editable = 'off';
            app.ErrorReport.FontColor = [1 0 0];
            app.ErrorReport.Position = [53 15 627 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Carriers

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end