function[fitObject] = fitBoseEinstein(x, y)
%in units of the band gap (meV typically)
aGuess = max(y) + 0.061;
bGuess = 0.061;
thetaGuess = 240;
fullGuess = [aGuess, bGuess, thetaGuess];

ft = fittype('boseEinstein(x, a, b,theta)');
[fitObject,gof,output] = fit(x,y,ft, 'Startpoint', fullGuess);
end