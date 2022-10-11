%% Change varshni parameters to bose-einstein parameters
E0 = 496.2; %meV
alpha = 2.8e-4*1000; %meV/K
beta = 300; %Kelvin
temps = linspace(4,300, 15)';

varshni = E0 - alpha*temps.^2./(temps + beta);

%plot(temps, varshni


%% 
fobject = fitBoseEinstein(temps, varshni);
plot(temps, varshni,'DisplayName', 'Varshni');
hold on
plot(temps, fobject(temps),'DisplayName', "Bose")
legend()
%%
e1 = 18;
e2 = 0.25;
n = 1/sqrt(2)* sqrt((e1 + sqrt((e1^2 + e2^2))));
k = 1/sqrt(2) * sqrt(-e1 + sqrt(e1^2 +e2^2));
lambda = 1.40254*1e-4;

alpha = 4*pi*k/lambda;