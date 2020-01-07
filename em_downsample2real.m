function [td,xd,yd] = em_downsample2real(t,x,y,sr)
%em_downsample2real  - downsamples eye position data to orignal (real) eyetracker sampling rate
%
% USAGE:
%	[t,x,y] = em_downsample2real(t,x,y,par.Downsample2Real);
%
% INPUTS:
%		t		- time axis (seconds)
%		x		- hor position (deg)
%		y		- ver position (deg)
%		sr		- real sampling rate
% OUTPUTS:
%		td,xd,yd	- downsampled data
%
% REQUIRES:	Neuroelf, Igtools
%
% See also em_smooth, em_filter
%
%
% Author(s):	I.Kagan, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2015-10-01:	Created function (Igor Kagan)
% 2020-01-07:   Changed removing consecutive identical segments to 2 x t_segment duration (from t_segment duration)
%               to account for variations in eyetracker sampling rate (IK)
% $Revision: 1.1 $  $Date: 2020-01-07 17:43:24 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

td = t;
xd = x;
yd = y;

t = t(:)';
x = x(:)';
y = y(:)';



t_segment = 2 * 1/sr; % s 
% we want to remove all consecutive identical segments of up to 2 x t_segment duration
% longer segments (e.g. due to signal dropout/blink saturation) will not be removed


dx = diff(x);
dy = diff(y);
repeats_idx = find(dx==0 & dy==0)+1;
if ~isempty(repeats_idx),
	repeated_segments_start_idx = [repeats_idx(1) repeats_idx(find(diff(repeats_idx)>1)+1)];
	repeated_segments_end_idx = [repeats_idx(find(diff(repeats_idx)>1)) repeats_idx(end)];

	segments2remove_idx = find((t(repeated_segments_end_idx) - t(repeated_segments_start_idx))<t_segment);

	idx2remove = cell2mat(arrayfun(@colon,repeated_segments_start_idx(segments2remove_idx),repeated_segments_end_idx(segments2remove_idx),'UniformOutput',false));

	td(idx2remove) = [];
	xd(idx2remove) = [];
	yd(idx2remove) = [];

end
if 0, % FOR DEBUG
	figure;
	plot(t,x,'g'); hold on; plot(t,x,'g.');
	plot(t,y,'m'); hold on; plot(t,y,'m.');
	plot(t(repeats_idx),x(repeats_idx),'rx');
	ig_add_multiple_vertical_lines(t(repeated_segments_start_idx),'Color','r');
	ig_add_multiple_vertical_lines(t(repeated_segments_end_idx),'Color','k');
	plot(t(idx2remove),x(idx2remove),'ro');
	plot(td,xd,'go');
end





