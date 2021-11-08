%% Plot the results
%%%%%% Here you need to double click on the sample %%%
%%%%%% "L##-###_FittedResults.mat" file and hit run for the script %%%
%%%%%% It will assume 9 temperature measurements were performed
%%%%%% So it will create a 3-by-3 matrix of TRPL decay data
%%%%%% A new figure is created for each excitation


n = 3; %%%Subplots matrix number of rows
m = 5; %%%Subplots matrix number of columns
test = ['b-', 'g-'];

 
% Create an array of figures 
figures = gobjects(length(Excitations) +1 ,1); 


%set this as the last figure
excitationsFigure = figures(end);


for exciteind = 1: length(Excitations)

 

figurename = string(SampleName) + '_AngleSetting_' + string(Excitations(exciteind));
newfig = figure('Name',figurename, 'NumberTitle', 'off');

figures(exciteind) = newfig;

    
    for tempind = 1:length(Temperatures) 
    %%make a subplot for each temperature
    subplot(n, m, tempind)
    hold on
    
    %%plot the Decay
    semilogy(Time_Axis, Decay(:,tempind, exciteind),'b.') 
    
    %%highlight the fit window in the decay
    plot(Time_Axis, Decay_WithinFitWindow(:,tempind, exciteind),'r.') 
    
    %%Plot the fitted decay
    plot(Time_Axis, FittedDecays(:,tempind, exciteind),'g-','LineWidth',1.5)
    
    %Denote the temperature and the characteristic slope
    title("T = " + string(Temperatures(tempind)) + " K, \tau = " + string(t0_vTemp(1,tempind, exciteind)) + "\mus")
    hold off
    set(gca,'YScale','log');
    xlim([0,1.5])
    end
     
%saveas(gcf, char(figurename), 'fig');
%saveas(gcf, char(figurename), 'png');

end

%% Plot the lifetimes

hold on
ColorWheel = [ 'bo-' ; 'go-' ; 'ro-' ; 'ko-' ; 'co-' ; 'yo-' ];

figurename = string(SampleName) + '_FinalResults';
newfig = figure('Name',figurename, 'NumberTitle', 'off');

figures(end) = newfig;
hold on
for exciteind = 1:length(Excitations)
   name = num2str(Excitations(exciteind));
   plot(Temperatures, t0_vTemp(1,:,exciteind), ColorWheel(exciteind,:),'Displayname',name ) 
 
    
end
legend
title(SampleName)
xlabel('Temperature (K)')
ylabel('Characteristic slope (\mus)')
%% Plot a single lifetime

%%% Double click on a FittedResults mat file and hit run!

%figure('Name','Comparison')

% the data name to the legend
dataName = 'Low Pass filter';

temperature = 77;  %Select the temperature that you want!
tempIndex = Temperatures == temperature; %this will get that temperature

subplot(2, 1, 1,'LineWidth', 2);
excitation = 39;   %Select the excitation that you want!
excitIndex = Excitations == excitation; %this will get that excitation
hold off
semilogy(Time_Axis, Decay(:,tempIndex, excitIndex), 'b--','DisplayName', dataName,'Linewidth', 1.3)
hold on
plot(Time_Axis, Decay_WithinFitWindow(:, tempIndex, excitIndex),'r.', 'MarkerSize', 8)
plot(Time_Axis, FittedDecays(:, tempIndex, excitIndex),'g-', 'LineWidth',2)








