% This is script changes the Data Fitting Status to
% WRITING

% Use just prior to data write to text file step, since it takes quite a
% while and we want users to know data is being written

% Tool status to WRITING, blue font
app.ToolStatus.Value = 'WRITING';
app.ToolStatus.FontColor = 'b'; 

% Pause ensures Tool Status updates before writing starts
pause(0.5)
