function out = em_saccade_blink_detection(t,x,y,varargin)
%em_saccade_blink_detection  - detects saccades and blinks from eye position data
%
% USAGE:
%	out = em_saccade_blink_detection(t,x,y,varargin);
% 	out = em_saccade_blink_detection(trial(k).tSample_from_time_start,trial(k).x_eye,trial(k).y_eye,'ma1_em_settings_monkey_220Hz');
%	out = em_saccade_blink_detection(t,x,y,'OpenFigure',true,'Plot',true)
%
% INPUTS:
%		t		- time axis (seconds)
%		x		- hor position (deg)
%		y		- ver position (deg)
%		varargin	- can be struct with (fieldname - value) pairs OR multiple (field)name - value pairs
%				- OR full path to m-file (not function!) containing settings variable as shown in custom_settings/em_custom_settings_example
%				- see default parameters below for available paramters 
% OUTPUTS:
%		out		- structure containing detected saccade parameters
%
% REQUIRES:	Neuroelf, Igtools, em
%
% See also em_smooth, em_filter
%
%
% Author(s):	I.Kagan, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2015-08-07:	Created function (Igor Kagan)
% 2015-10-01:	Added option to downsample the traces to "real" (e.g. 220Hz, or 60Hz) sampling rate before interpolation (IK): em_downsample2real.m
% ...
% $Revision: 1.0 $  $Date: 2015-10-01 17:43:24 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%



% default parameters
defpar = { ...
	'SampleRate',		'double',	'nonempty',	1000; ...		% Hz
	'Downsample2Real',	'double',	'nonempty',	0; ...			% Hz (actual eyetracker sampling rate), if 0 - do not downsample 
	'SacOnsetVelThr',	'double',	'nonempty',	50; ...			% deg/s
	'SacOffsetVelThr',	'double',	'nonempty',	20; ...			% deg/s
	'MinSacDuration',	'double',	'nonempty',	0.02; ...		% s
	'MaxSacDuration',	'double',	'nonempty',	Inf; ...		% s
	'MinSacAmplitude',	'double',	'nonempty',	0.25; ...		% deg
	'MaxSacAmplitude',	'double',	'nonempty',	Inf; ...		% deg
	'PosSmoothConvWin',	'char',		'nonempty',	'gausswin'; ...		% 'rectwin', 'gausswin', etc. see 'help window'
	'PosSmoothConvLen',	'double',	'nonempty',	0; ...			% s, length of conv kernel, set to 0 of no smooting
	'PosFilterCutoff',	'double',	'nonempty',	30; ...			% Hz, set to 0 if no filter	
	'VelSmoothConvWin',	'char',		'nonempty',	'gausswin'; ...		% 'rectwin', 'gausswin', etc. see 'help window'
	'VelSmoothConvLen',	'double',	'nonempty',	0; ...			% s, length of conv kernel, set to 0 of no smooting
	'VelFilterCutoff',	'double',	'nonempty',	30; ...			% Hz, set to 0 if no filter
	'VelAdaptiveThr',	'logical',	'nonempty',	false; ...		% true or false
    'MinFixDurAfterSac','double',	'nonempty',	0; ...          % s
	
	'Plot',		'logical',	'nonempty',	false; ...			% true or false
	'OpenFigure',	'logical',	'nonempty',	false; ...			% true or false
	'Verbose',	'logical',	'nonempty',	false; ...			% true or false
	
	'xposcol',	'double',	'nonempty',	[0 0.4980 0]; ...
	'yposcol',	'double',	'nonempty',	[0.4784 0.0627 0.8941]; ...
	'sxposcol',	'double',	'nonempty',	[0 1 0]; ...
	'syposcol',	'double',	'nonempty',	[1 0 1]; ...	
	'velcol',	'double',	'nonempty',	[0 0 0.5]; ...
	'svelcol',	'double',	'nonempty',	[0 0 1]; ...
	'sac_oo_col',	'double',	'nonempty',	[0.7294    0.8314    0.9569]; ...
	};

if nargin > 3, % specified parameters
	if ~isstruct(varargin{1}),
		if length(varargin)==1, % specified .m file with settings, such as em_custom_settings_example.m
			run(varargin{1});
			par = struct(settings{:});
		else % parameter-value pair(s)
			par = struct(varargin{:});
		end
	else
		par = varargin{1};
	end
	par = checkstruct(par, defpar);
else
	par = checkstruct(struct, defpar);
end

t = double(t);
x = double(x);
y = double(y);

original_SI = mode(diff(t)); % original sampling interval

if par.Downsample2Real,
	if par.Verbose,
		disp(sprintf('Downsampling from %.5f Hz to %.5f Hz',1/original_SI,par.Downsample2Real));
	end
	% the idea is to take the first sample of "artifically" (e.g. due to monkeypsych) upsampled data
	[t,x,y] = em_downsample2real(t,x,y,par.Downsample2Real);
end

if length(unique(round2(diff(t),0.0001))) > 1 || abs(original_SI-1/par.SampleRate)>0.001, % resample
	if par.Verbose,
		disp(sprintf('Interpolating from %.5f Hz to %.5f Hz',1/original_SI,par.SampleRate));
	end
	ti = t(1) : 1/par.SampleRate : t(end);
	xi = interp1(t, x, ti, 'linear');
	yi = interp1(t, y, ti, 'linear');
else
	% data in rows
	ti = t(:)';
	xi = x(:)';
	yi = y(:)';
end


% Position
if par.PosSmoothConvLen
	if par.Verbose,
		disp(sprintf('Smoothing position with %s of %d samples',par.PosSmoothConvWin,par.SampleRate*par.PosSmoothConvLen));
	end
	xis = em_smooth(xi,par.PosSmoothConvWin,par.SampleRate*par.PosSmoothConvLen);
	yis = em_smooth(yi,par.PosSmoothConvWin,par.SampleRate*par.PosSmoothConvLen);
end

if par.PosFilterCutoff
	if par.Verbose,
		disp(sprintf('Filtering position at %d Hz using %d Hz sampling rate ',par.PosFilterCutoff,par.SampleRate));
	end
	
	if exist('xis','var');
		xis = em_filter(xis,par.PosFilterCutoff,par.SampleRate);
		yis = em_filter(yis,par.PosFilterCutoff,par.SampleRate);
	else
		xis = em_filter(xi,par.PosFilterCutoff,par.SampleRate);
		yis = em_filter(yi,par.PosFilterCutoff,par.SampleRate);
	end
end

% Velocity
if exist('xis','var'),
	v = [0 sqrt( (diff(xis).*par.SampleRate).^2+(diff(yis).*par.SampleRate).^2 )];
else
	v = [0 sqrt( (diff(xi).*par.SampleRate).^2+(diff(yi).*par.SampleRate).^2 )];
end

if par.VelSmoothConvLen
	if par.Verbose,
		disp(sprintf('Smoothing velocity with %s of %d samples',par.VelSmoothConvWin,par.SampleRate*par.VelSmoothConvLen));
	end
	
	vs = em_smooth(v,par.VelSmoothConvWin,par.SampleRate*par.VelSmoothConvLen);
else
	vs = v;
end

if par.VelFilterCutoff
	if par.Verbose,
		disp(sprintf('Filtering velocity at %d Hz using %d Hz sampling rate ',par.VelFilterCutoff,par.SampleRate));
	end
	vs = em_filter(vs,par.VelFilterCutoff,par.SampleRate);	
end

% Velocity threshold crossings
[on_i off_i] = ig_find_threshold_crossings(vs,[par.SacOnsetVelThr par.SacOffsetVelThr]);

% Duration
dur = ti(off_i) - ti(on_i);
idx_dur_valid = find(dur > par.MinSacDuration & dur < par.MaxSacDuration);

% Amplitude
amp = sqrt( (xi(on_i)-xi(off_i)).^2 + (yi(on_i)-yi(off_i)).^2);
idx_amp_valid = find(amp > par.MinSacAmplitude & amp < par.MaxSacAmplitude);

valid_idx = intersect(idx_dur_valid,idx_amp_valid);

sac_onsets	= ti(on_i(valid_idx));
sac_offsets	= ti(off_i(valid_idx));

% Minimal fixation period after saccade
if par.MinFixDurAfterSac,
    FixDur = sac_onsets(2:end) - sac_offsets(1:end-1);
    ind_fix_dur_valid = [find(FixDur > par.MinFixDurAfterSac) length(valid_idx)];
    valid_idx = valid_idx(ind_fix_dur_valid);
    sac_onsets = sac_onsets(ind_fix_dur_valid);
    sac_offsets = sac_offsets(ind_fix_dur_valid);
end

out.par = par;
out.sac_onsets = sac_onsets;
out.sac_offsets = sac_offsets;
out.sac_amp = amp(valid_idx);
out.sac_dur = dur(valid_idx);

for k = 1:length(valid_idx)
	out.sac_max_vel(k)	= max(vs(on_i(valid_idx(k)):off_i(valid_idx(k))));
	out.sac_mean_vel(k)	= mean(vs(on_i(valid_idx(k)):off_i(valid_idx(k))));
end


if par.Plot,
	
	if par.OpenFigure,
		figure('Name','em_saccade_blink_detection');

		
		subplot(2,1,1); hold on;
		title(linewrap(ig_struct2string(par,{'SampleRate' 'Downsample2Real' 'SacOnsetVelThr' 'SacOffsetVelThr' 'MinSacDuration' 'MaxSacDuration' ...
			'MinSacAmplitude' 'MaxSacAmplitude' 'PosFilterCutoff' 'VelFilterCutoff'}),128),'FontSize',8);
		
		plot(ti,xi,'k-','Color',par.xposcol);
		plot(ti,yi,'k-','Color',par.yposcol);
		ylim = get(gca,'Ylim');
		for k=1:length(sac_onsets),
			patch([sac_onsets(k) sac_offsets(k) sac_offsets(k) sac_onsets(k)],[ylim(2) ylim(2) ylim(1) ylim(1)],par.sac_oo_col);
		end
		plot(t,x,'k.','Color',par.xposcol);
		plot(ti,xi,'k-','Color',par.xposcol);
		plot(t,y,'k.','Color',par.yposcol);
		plot(ti,yi,'k-','Color',par.yposcol);
		if exist('xis','var'),
			plot(ti,xis,'k-','Color',par.sxposcol);
			plot(ti,yis,'k-','Color',par.syposcol);
		end
		ylabel('Eye position');

		subplot(2,1,2); hold on;
		plot(ti,vs,'k-','Color',par.velcol);
		ylim = get(gca,'Ylim');
		for k=1:length(sac_onsets),
			patch([sac_onsets(k) sac_offsets(k) sac_offsets(k) sac_onsets(k)],[ylim(2) ylim(2) ylim(1) ylim(1)],par.sac_oo_col);
		end

		% for debug
		plot(ti,v,'k.','Color',par.velcol);
		plot(ti,v,'k-','Color',par.velcol);

		if exist('vs','var'),
			plot(ti,vs,'k.','Color',par.svelcol);
			plot(ti,vs,'k-','Color',par.svelcol);
		else
			plot(ti,v,'k.','Color',par.velcol);
			plot(ti,v,'k-','Color',par.velcol);
		end

		line([ti(1) ti(end)],[par.SacOnsetVelThr par.SacOnsetVelThr],'Color',[1 0 0]);
		line([ti(1) ti(end)],[par.SacOffsetVelThr par.SacOffsetVelThr],'Color',[ 0.8471    0.1608         0]);
		ylabel('Eye velocity');

		ig_set_all_axes('Xlim',[ti(1) ti(end)]);
	
	else % just plot saccade detection results on the current axis
		
		ig_add_multiple_vertical_lines(sac_onsets,'Color',[1 0 0],'LineStyle',':');
		ig_add_multiple_vertical_lines(sac_offsets,'Color',[0.6000    0.2000         0],'LineStyle',':');
		
		
	end
		
end
