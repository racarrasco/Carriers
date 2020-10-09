function[totalLifetime] = calculateTotalLifetimeLink(...
    x, a, b, c, d, e, f, g,...
    h, k, l, m1, m2, m3, n1,...
%{
n2,
    %}
n3, o1, o2,o3)
% Calculate all of the lifetime components and add them together. The rates
% are added then the reciprocal of the total rate is the total lifetime
% This is for the fitting routine where the single output is the total
% lifetime output, SO THIS FILE WILL BE HARD TO READ, IN ORDER TO
% UNDERSTAND THE CODE BETTER (AND QUICKER) SEE THE FILE 
% calculateLifetimes.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Input Parameters %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x = Temperature that will be used to calculate the total lifetime
% a =  Temperature lengths of each sample
% b = doping type n-type or p-type
% c = egTempDep the temperature dependent band gap
% d = Ev the valence edge energy (taken as the negative band gap)
% e = Ec the conduction edge energy (taken as the 0 reference) 
% f = einf the static dielectric constant below the absorption edge
% g = phi the photon recycling factor



% FITTING PARAMETERS
% h = mestar 
% k = mhstar
% l = Auger overlap parameter  OPTIONAL
% m = doping density in units of 1e-15cm^(-3)
% n = defect level in reference to the conduction band
% o = defect density in units of inverse meters
totalLifetime = NaN*ones(size(x));

xone = x(1:a(1));
c1 = c(1:a(1));
d1 = d(1:a(1));


if (size(a,1) > 1)

xtwo = x(a(1) + 1: a(1)+a(2));
c2 = c(a(1) + 1:a(1) + a(2));
d2 = d(a(1) + 1:a(1) + a(2));



lastIndex = sum(a(1:2));
xthree = x(lastIndex + 1: lastIndex+a(3));
c3 = c(lastIndex + 1: lastIndex+a(3));
d3 = d(lastIndex + 1: lastIndex+a(3));
totalLifetime(a(1) + 1: a(1)+a(2)) = calculateTotalLifetime3(xtwo, b(2), c2', d2', e, f,g(2),...
    h, k, l,m2,n1,o2);

totalLifetime(lastIndex + 1: lastIndex+a(3)) = calculateTotalLifetime4(xthree, b(3), c3', d3', e, f,g(3),...
    h, k, l,m3,n1,o1, n3, o3);

end
%You will need to change this!
totalLifetime(1:a(1)) = calculateTotalLifetime3(xone, b(1), c1', d1', e, f,g(1),...
    h, k, l,m1,n1,o1);






end
