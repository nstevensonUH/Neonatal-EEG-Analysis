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

load  kirsi_BD_filters_256 % HP 6th order Elliptical fc = 10Hz. LP 1st order Butterworth fc = 0.5Hz

% Initialise Parameters
fs1 = 256; fs2 = 16;
dat = resample(eeg_data, fs1, Fs); % NB - must be at 256 as NLEO is different when when fs changes so parameters above will need to change is fs is changed.

dat = filter(Num_LP, Den_LP, dat);   % Filter data between 0.5-10Hz
dat = filter(Num_HP, Den_HP, dat);

snleo = nlin_energy(dat, win*fs1);      % Smoothed Absolute Nonlinear Energy Operator
s1 = resample(snleo, fs2, fs1); s2 = s1;
epl = fs2*wc; 
for ii = epl+1:length(s1);                      % Baseline Correction
    r1 = ii-epl; r2 = ii-1; 
    s2(ii) = s1(ii)-min(s1(r1:r2));
   % s3(ii) = min(s1(r1:r2));
end

% M = 1000; N = length(s2);
% th = linspace(quantile(s2, 0.1), quantile(s2, 0.9), M);
% erf1 = zeros(1, length(th)-1);
% for ii = 1:length(th)-1;
%     dum = zeros(1, N);
%     dum(s2>th(ii)) = 1;
%     r1 = find(diff([0 dum 0]) == 1);
%     erf1(ii) = length(r1);   % This is the cost function which threshold selection is based
% end
% q1 = find([erf1 0]<[0 erf1]);  % indices where erf1 decreases
% % A few checks and if not satisfied merely select the threshold as the
% median of the data
% The aim is to select the start of the longest period of constant (stable)
% ERF which is just burst number slected by a given threshold
% if isempty(q1)==1
%     th1 = th(M/2);
% else
%     [q2, q3]= sort(diff(q1), 'descend');
%     q4 = find(erf1(q1(q3))>1, 1);
%     if isempty(q4)==0; 
%         th1 = th(q1(q3(q4))); 
%     else
%         th1 = th(M/2);
%     end
% end
% dum = zeros(1, N);
% dum(snleo>th1) = 1;
% t = linspace(0,8, length(dat));
% plot(t, snleo); hold on;
% plot(t, th1*dum)


bm = zeros(1, length(s1));                  % Generate Burst Mask
bm(s2>th) = 1;
% r1 = find(diff([0 bm 0])==1);
% r2 = find(diff([0 bm 0])==-1);
% for ii = 1:length(qq)                               % Remove short duration bursts
%     bm(r1(qq(ii)):r2(qq(ii))-1)=0;
% end
bm = process_ba_1ch(bm, fs2, 1, 2, 0);
%t = 0:1/fs2:(length(bm)-1)/fs2;

ba = zeros(1,length(bm)*4);
for ii = 1:length(bm);
    q1 = 1+(ii-1)*4; q2 = q1+3;
    ba(q1:q2) = bm(ii).*ones(1,4);
end

if length(ba)~=length(eeg_data);
    ba = [ba zeros(1,length(eeg_data)-length(ba))];
end

end

