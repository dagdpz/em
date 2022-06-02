% em_custom_settings_monkey_220Hz
% settings file for use for em_saccade_blink_detection
% see defpar in em_saccade_blink_detection.m for more options
% should contain variable **settings** as below!

settings = {...
	
	'SampleRate',		1000, ...		% Hz
	'Downsample2Real',	220, ...		% Hz (actual eyetracker sampling rate), if 0 - do not downsample 
	'SacOnsetVelThr',	30, ...			% deg/s
	'SacOffsetVelThr',	10, ...			% deg/s
	'MinSacDuration',	0.03, ...		% s
	'MaxSacDuration',	Inf, ...		% s
	'MinSacAmplitude',	0.5, ...         % deg
	'MaxSacAmplitude',	Inf, ...		% deg
	'PosSmoothConvWin',	'gausswin', ...	% 'rectwin', 'gausswin', etc. see 'help window'
	'PosSmoothConvLen',	0, ...			% s, length of conv kernel, set to 0 of no smooting
	'PosFilterCutoff',	100, ...			% Hz, set to 0 if no filter	
	'VelSmoothConvWin',	'gausswin', ...	% 'rectwin', 'gausswin', etc. see 'help window'
	'VelSmoothConvLen',	0, ...			% s, length of conv kernel, set to 0 of no smooting
	'VelFilterCutoff',	100, ...			% Hz, set to 0 if no filter
	'VelAdaptiveThr',	false, ...		% true or false
    'MinFixDurAfterSac',0, ...          % s
	'Plot',             true, ...       % true or false (if true, results of saccade detection will be plotted)
	'OpenFigure',		true, ...		% true or false (if true, new figure showing saccade detection will be open)
	
	};