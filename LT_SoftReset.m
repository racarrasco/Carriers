% This is the soft reset script which serves to partially reset the LT user interface

% The soft reset will only set the Tool Status to INACTIVE, and will leave
% all other options/settings alone.  This is intended to be used after an
% error is identified, so the user can see the settings prior to
% encountering the error

% Tool status to INACTIVE, red font
app.ToolStatus.Value = 'INACTIVE';
app.ToolStatus.FontColor = 'r'; 
