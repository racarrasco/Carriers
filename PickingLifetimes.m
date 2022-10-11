%% pick out specific temperatures

tempsInterest = 70:10:300;

lifetimesInterest = nan.*ones(length(tempsInterest),1);

for index = 1:length(tempsInterest)
    tempIndex = tempsInterest(index) == Tprobe;
    lifetimesInterest(index) = totalLifetime(tempIndex);
    
end


plot(Tprobe, totalLifetime)
hold on
plot(tempsInterest, lifetimesInterest,'o')