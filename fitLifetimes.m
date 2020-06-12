function[fitObject, gof, output] = ...
 fitLifetimes(xin,yin, fitAugerOverlap,fixSRH1, type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
  einf, Nc, Nv, ni,G, phi, f1f2, dopingDensity, defectLevel, defectDensity,...
  activation1, defectLevel2, defectDensity2, activation2)


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




if(nargin == 21) 
    if(fitAugerOverlap) 
            % CASE 1
            % Fit auger overlap, doping, defect level and defect density
            ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r,s)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n']);
            guess = [f1f2, dopingDensity, defectLevel, defectDensity, activation1]; 
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi},...
            'Lower', [0.01 , 0.00001, -1, 0, 0],...
            'Upper', [0.3 , 100000, 0, inf, 1], ...
            'Startpoint', guess);  

     
    else
        % CASE 2
        % Only fit 3 parameters which are doping Density, defect Level, and
        % defect density, NOT bloch overlap
        
        ft = fittype('calculateTotalLifetime(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s)',...
        'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o']);
        % The guess is doping density, defect Level and defect density
        guess = [dopingDensity, defectLevel, defectDensity, activation1];
     
        [fitObject, gof, output] = fit(xin',yin',ft,...
        'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
        einf, Nc, Nv, ni, G, phi, f1f2},...
        'Lower', [0.000001, -1, 0,0],...
        'Upper', [100000, 0, inf, 1], ...
        'Startpoint', guess);

    end 
    
    
    
    
    
    
    
% We have more input arguments for fitting, which means that we have intentional
% doping
else
    if (fitAugerOverlap)
        if(fixSRH1) 
        %FIX SRH1 parameters, fit SRH2 parameters, Bloch overlap parameter,
        %doping (Fit 4 parameters)
            ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t, u, v)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'q'; 'r'; 's']);
        
        %Fit parameter o(Bloch overlap), p (doping density), s(defectLevel2), t(defectDensity2) 
             guess = [f1f2, dopingDensity, defectLevel2, defectDensity2, activation2];
        
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi,defectLevel, defectDensity, activation1},...
            'Lower', [0.01 , 0.00001, -1, 0],...
            'Upper', [0.3 , 100000, 0, inf], ...
            'Startpoint', guess);    
        else
        % Fit SRH1, SRH2, Bloch overlap parameter, doping (8 parameters)
            ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t, u, v)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n']);
        
        %Fit parameter o(Bloch overlap) p(doping density) q(defectLevel1) r(defectDensity1)
        %s (activation1) t(defect level2) u(defect density2) v
        %(activation2)
            guess = [f1f2, dopingDensity, defectLevel, defectDensity, activation1,...
                defectLevel2, defectDensity2, activation2];
            
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi},...
            'Lower', [0.01 , 0.00001, -1, 0, -1, 0],...
            'Upper', [0.3 , 100000, 0, inf, 0, inf], ...
            'Startpoint', guess);  
            
        
        
        
        
            
        end
            
    
    else 
        if(fixSRH1)
             %FIX SRH1 and Bloch overlap. 
            
             %FIX SRH1 parameters and Bloch overlap, fit SRH2 parameters,
        %doping (Fit 4 parameters)
            ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t, u, v)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n'; 'o'; 'q'; 'r';'s']);
        
        %Fit parameter  p (doping density), t(defectLevel2), u(defectDensity2) v(activation2) 
             guess = [dopingDensity, defectLevel2, defectDensity2, activation2];
        
             
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi,f1f2, defectLevel, defectDensity,...
            activation1},...
            'Lower', [0.00001, -1, 0,0],...
            'Upper', [10000, 0, inf,1], ...
            'Startpoint', guess);
            
        else
            
            
            % Fit SRH1, SRH2, doping, FIX Bloch overlap parameter (5 parameters)
            ft = fittype('calculateTotalLifetime2(x, a, b, c, d, e, f, g, h, k, l, m, n, o, p, q, r, s, t)',...
            'problem',['a'; 'b'; 'c'; 'd'; 'e'; 'f'; 'g'; 'h'; 'k'; 'l'; 'm'; 'n';'o']);
        
        %Fit parameter  p(doping density) q(defectLevel1) r(defectDensity1)
        % s(defect level2) t(defect density2)
            guess = [ dopingDensity, defectLevel, defectDensity,...
                activation1,defectLevel2, defectDensity2, activation2];
            
            
            [fitObject, gof, output] = fit(xin',yin',ft,...
            'problem',{type, meStar, mhStar, eg, valenceEdge, conductionEdge,...
            einf, Nc, Nv, ni, G, phi, f1f2},...
            'Lower', [ 0.00001, -1, 0,0, -1, 0,0],...
            'Upper', [ 100000, 0, inf,1, 0, inf,1], ...
            'Startpoint', guess);  
               
       
        end
        
        
        
    end
        
end








