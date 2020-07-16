function[fitObject, gof, output] = ...
 fitLifetimes2(xin,yin, fitAugerOverlap,fixSRH1, type, eg, valenceEdge, conductionEdge,...
  einf, phi, meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity,...
  defectLevel2, defectDensity2)


% Calculate all of the lifetime components and add them together the rates
% are added then the reciprocal of the total rate is the total lifetime
% SEE THE calculateLifetimes.m file for meaning to the input parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%yin   =  the data we will be performing on the fitting

%These will be the fit parameters, and MATLAB arranges the fit parameters
%alphabetically, not in the order in which they appear.


% There are a lot of forks due to the multiple checkboxes
% Case 1: fit Bloch overlap, doping, defect1 level and defect1 density,
% there is no secondary SRH component

% Case 2: nargin == 20, but Bloch checkbox is false
% fit doping, defect1 level, and defect1 density but do not fit Bloch overlap 

% Give a weighting to each measurement and just use that value as the
% weight
% 1 us error is 1us,   and 400 ns error is 400 ns etc.
%w = [];
w = yin';




if(nargin == 16) 
    if(fitAugerOverlap) 
            % CASE 1
            % Fit Effective masses auger overlap, doping,
            % defect level and defect density,
            ft = fittype('calculateTotalLifetime3(x, a, b, c, d, e, f, g, h, k, l, m, n)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f']);
            guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity]; 
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge, einf, phi},...
            'Lower', [0.001, 0.001, 0.01, 0.00001, -1, 0],...
            'Upper', [ 2,        inf,  0.3, 1000, 0, inf], ...
            'Startpoint', guess,...
            'Weights', w);  

     
    else
        % CASE 2
        % Only fit 5 parameters which are effective masses (2), doping Density,
        % defect Level, and defect density,NOT bloch overlap
        
        ft = fittype('calculateTotalLifetime3(x, a, b, c, d, e, f, g, h, k, l, m, n)',...
        'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k';]);
        % The guess is doping density, defect Level and defect density
        guess = [meStar, mhStar,dopingDensity, defectLevel, defectDensity];
     
        [fitObject, gof, output] = fit(xin',yin',ft,...
        'problem',{type, eg, valenceEdge, conductionEdge,...
        einf, phi, f1f2},...
            'Lower', [0.001, 0.001, 0.00001, -1, 0],...
            'Upper', [ 2,        inf, 1000, 0, inf], ...
        'Startpoint', guess,...
        'Weights', w);

    end 
    
    
    
    
    
    
    
% We have more input arguments for fitting, which means that we have intentional
% doping
else
    if (fitAugerOverlap)
        if(fixSRH1) 
        %FIX SRH1 parameters, fit effective masses
        %SRH2 parameters, Bloch overlap parameter,
        %doping (Fit 6 parameters)
            ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'm'; 'n']);
        
        %Fit parameter effective masses (g,h) k(Bloch overlap),
        %l (doping density), o(defectLevel2), p(defectDensity2) 
           
             guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel2, defectDensity2];
        
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge,...
            einf, phi,defectLevel, defectDensity},...
            'Lower', [0.001, 0.001, 0.01, 0.00001, -1, 0],...
            'Upper', [2,         inf,  0.3, 100000, 0, inf], ...
            'Startpoint', guess,...
            'Weights', w);    
        else
        % Fit SRH1, SRH2, Bloch overlap parameter, doping and effective maases (8 parameters)
            ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f']);
        
        %Fit parameter effective mass (g,h) k(Bloch overlap) l(doping density)
        % m(defectLevel1) n(defectDensity1) o(defect level2) p(defect density2)
            guess = [meStar, mhStar, f1f2, dopingDensity, defectLevel, defectDensity,...
                defectLevel2, defectDensity2];
            
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge,...
            einf, phi},...
            'Lower', [0.001, 0.001,0.01 , 0.00001, -1, 0, -1, 0],...
            'Upper', [2,         inf,0.3 , 100000, 0, inf, 0, inf], ...
            'Startpoint', guess,...
            'Weights', w);  
            
        
        
        
        
            
        end
            
    
    else 
        if(fixSRH1)
             %FIX SRH1 and Bloch overlap. 
            
             %FIX SRH1 parameters and Bloch overlap, fit SRH2 parameters,
        %doping and effective masses (Fit 5 parameters)
            ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k'; 'm'; 'n']);
        
        %Fit parameter  p (doping density), s(defectLevel2), t(defectDensity2) 
             guess = [meStar, mhStar, dopingDensity, defectLevel2, defectDensity2];
        
             
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge,...
            einf,  phi,f1f2, defectLevel, defectDensity},...
            'Lower', [0.001, 0.001,0.00001, -1, 0],...
            'Upper', [2,         inf,100000, 0, inf], ...
            'Startpoint', guess,...
            'Weights', w);
            
        else
            
            
            % Fit effective mass, SRH1, SRH2, doping, FIX Bloxh overlap parameter (7 parameters)
            ft = fittype('calculateTotalLifetime4(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'k']);
       
        %Fit parameter gh (effective masses) l(doping density) m(defectLevel1) n(defectDensity1)
        % o(defect level2) (defect density2)
            guess = [meStar, mhStar, dopingDensity, defectLevel, defectDensity,...
                defectLevel2, defectDensity2];
            
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, eg, valenceEdge, conductionEdge,...
            einf, phi, f1f2},...
            'Lower', [0.001, 0.001, 0.00001, -1, 0, -1, 0],...
            'Upper', [2,         inf, 100000, 0, inf, 0, inf], ...
            'Startpoint', guess,...
            'Weights', w);  
               
    
        end
        
        
        
    end
        
end