function ba = detector_per_channel_palmu(eeg_data, Fs, th1)
% Kirsi's EEG 'burst' detector. Code based on the paper - Palmu et al,
% Physiol Meas 31:85-93, 2010.
%
% [t, bm] = detector_per_channel_palmu(eeg_data, Fs)
%
% Inputs:- a single channel of EEG (vector 1xN)
%              - sampling frequency of EEG
%
% Outputs - time vector (output is sampled at 16Hz - vector 1xM)
%                 - binary burst mask (1 - burst, 0 - not burst - vector 1xM)
%
% Nathan Stevenson
% August 2016
% University of Helsinki

% Initialise Parameters (feel free to tinker)
win = 1.5;
wc = 60;
md = 1;
th = 1.2/15*th1-0.1;%1.5;
if th<0.3; th = 0.3; end
if th>2.5; th = 2.5; end
load kirsi_BD_filters_256 % HP 6th order Elliptical fc = 10Hz. LP 1st order Butterworth fc = 0.5Hz
% Initialise Parameters
fs1 = 256; fs2 = 16;
dat = resample(eeg_data, fs1, Fs); % NB - must be at 256 as NLEO is different when when fs changes so parameters above will need to change is fs is changed.
dat = filter(Num_LP, Den_LP, dat);   % Filter data between 0.5-10Hz
dat = filter(Num_HP, Den_HP, dat);
snleo = nlin_energy(dat, win*fs1);      % Smoothed Absolute Nonlinear Energy Operator
s1 = resample(snleo, fs2, fs1); s2 = s1;
epl = fs2*wc; 
for ii = epl+1:length(s1)                      % Baseline Correction
    r1 = ii-epl; r2 = ii-1; 
    s2(ii) = s1(ii)-min(s1(r1:r2));
end
bm = zeros(1, length(s1));                  % Generate Burst Mask
bm(s2>th) = 1;
bm = process_ba_1ch(bm, fs2, 1, 2, 0);
ba = zeros(1,length(bm)*(fs1/Fs));
for ii = 1:length(bm)
    q1 = 1+(ii-1)*4; q2 = q1+3;
    ba(q1:q2) = bm(ii).*ones(1,4);
end

if length(ba)~=length(eeg_data)
    ba = [ba zeros(1,length(eeg_data)-length(ba))];
end

end

