function s = em_filter(x,FilterCutoff,SampleRate)
%em_filter  - filters eye movement data
%
% USAGE:	
% usage example1;
%
% INPUTS:
%		x		- signal (position or velocity)
%		FilterCutoff	- Hz
%		SampleRate	- Hz
% OUTPUTS:
%		s		- filtered signal
%
% REQUIRES:	Signal Processing Toolbox (help signal)
%
% See also cheby2, em_smooth, em_saccade_blink_detection
%
%
% Author(s):	I.Kagan, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2015-08-07:	Created function (Igor Kagan)
% ...
% $Revision: 1.0 $  $Date: 2015-08-07 17:47:22 $

% ADDITIONAL INFO:
% SampleRate/2 corresponds to "1.0"
% in cheby2: 0.0 < Wst < 1.0, with 1.0 corresponding to half the sample rate.
% Thus, for example normalized stopband edge frequency cutoff at 30 Hz would be Wst = 0.06;
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%% 

Wst = FilterCutoff/(0.5*SampleRate);
[B, A] = cheby2(1,3, Wst);

s=filtfilt(B,A,x);

