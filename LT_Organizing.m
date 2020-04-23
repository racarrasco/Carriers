% This is script changes the Data Fitting Status to
% ORGANIZING

% Lets the user know that the tool is organizing a dataset

% Tool status to ORGANIZING, cyan font
app.ToolStatus.Value = 'ORGANIZING';
app.ToolStatus.FontColor = 'c'; 

% Pause ensures Tool Status updates before processes start
pause(0.5)
