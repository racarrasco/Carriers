function [ gamma ] = get_photon_recycling_factor( a , n , d )
% Gets the photon recycling factor

% Input Parameters
% a = Bandgap absorption coefficient                       (1/cm)
% n = Refractive index of cap layer at bandgap energy
% d = Active region thickness                              (nm)

% Output Parameters
% gamma = photon recycling factor

% The photon recycling factor is the sum of three double integrals (r1 + r2 + r3)
% divided by another integral (denom)

% First convert thickness d from nm into cm
d = d*1e-7;

% The variables being integrated are distance into active region x (rows) 
% and emission angle theta (columns)
dx = d/1000;
x = 0:dx:d;
x = x';
dtheta = (pi/2)/999;
theta = 0:dtheta:pi/2;

% Get reflection coefficient at surface
R = ((n-1)/(n+1))^2;

% Get critical angle
n_outside = 1; % Refractive index of vacuum
thetaC = asin( n_outside / n ); % radians

% The index of the critical angle within theta is
[ ~ , thetaC_index ] = min( abs( theta - thetaC ) );

% Repmat x and theta so that element-wise multiplications can be performed
L_x = length(x);
L_theta = length(theta);
x     = repmat( x , 1 , L_theta );
theta = repmat( theta , L_x , 1);

% r1 = rays emitted towards surface within escape cone
% 0 < theta < thetaC
% 0 < x < d
r1_integrand = ( 1 - exp(-a*x./cos(theta)) + R*(1-exp(-a*d./cos(theta))).*exp(-a*x./cos(theta)) ).*sin(theta).*cos(theta);

% r2 = rays emitted towards surface outside of escape cone
% thetaC < theta < pi/2
% 0 < x < d
r2_integrand = ( 1 - exp(-a*x./cos(theta)) + (1-exp(-a*d./cos(theta))).*exp(-a*x./cos(theta)) ).*sin(theta).*cos(theta);

% r3 = rays emitted towards substrate
% 0 < theta < pi/2
% 0 < x < d
r3_integrand = ( 1 - exp(-a*(d-x)./cos(theta)) ).*sin(theta).*cos(theta);

% Denominator of photon recycling factor expression
denom_integrand = sin(theta).*cos(theta);

% Integrate using trapz()
% display('r1')
r1_thetaIntegrated = dtheta*trapz( r1_integrand(:,1:thetaC_index) , 2 );
% size(r1_thetaIntegrated)
r1 = dx*trapz( r1_thetaIntegrated , 1 );
% size(r1)

% display('r2')
r2_thetaIntegrated = dtheta*trapz( r2_integrand(:,thetaC_index:L_theta) , 2 );
% size(r2_thetaIntegrated)
r2 = dx*trapz( r2_thetaIntegrated , 1 );
% size(r2)

% display('r3')
r3_thetaIntegrated = dtheta*trapz( r3_integrand , 2 );
% size(r3_thetaIntegrated)
r3 = dx*trapz( r3_thetaIntegrated , 1 );
% size(r3)

% Only need to integrate along one row of denom, since all rows are
% identical repmat rows
denom = 2*d*( dtheta*trapz( denom_integrand(1,:) , 2 ) );

% Photon recycling factor
gamma = ( r1 +  r2 + r3 )/denom;