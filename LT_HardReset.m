% This is the hard reset script which serves to reset the LT user interface

% The hard reset is a full reset of the interface

% Tool status to INACTIVE, red font
app.ToolStatus.Value = 'INACTIVE';
app.ToolStatus.FontColor = 'r'; 



% Error Report
app.ErrorReport.Value = ''; % Error report field - Clear

% User Specified Axis Limits Panel
% app.UserAxisLimits.Value = 0; % User axis limits checkbox - Uncheck

% Decay Plot - Clear the decay plot
title( app.UIAxes , 'Temperature Dependence' )
scatter( app.UIAxes , 1 , 1 , 'x' )  
xlim( app.UIAxes , 'Auto' )
ylim( app.UIAxes , 'Auto' ) 