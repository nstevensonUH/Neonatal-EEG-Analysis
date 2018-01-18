function  [bursts1, art, ndm, reeg, ch_left, ch_right] = estimate_features_preprocessing(dat, fs1, nref);
% This function loads and pre-processes the EEG data file in fname  
% Preprocessing stages include
%   1. filtering and resampling
%   2. rEEG estimation
%   3. burst/SAT detection
%   4. artefact detection (high amplitude and zero amplitude impedence checks)
%   5. sleep state detection
%   6. Electrode short detection
%
% INPUTS: fname - filename (string)
%                  ch_no - the number of channels recorded in the EEG file
%                  ch_h - a vector denoting which channels are on the left
%                  hemisphere (=1) and which are on the right (=2)
%
% OUTPUTS: bursts1 - matrix containing burst annotation (0 no burst, 1 burst)
%                      va - vector contatining sleep state annotation (0 quiet sleep, 1 active sleep)
%                      art - vector containing artefac annotation (0 no artefact, 1 artefact)
%                      ndm - matrix containg the EEG data for each channel
%                      reeg - matrix containing the rEEG data for each channel (will be shorter than ndm as rEEG is sampled every 2s)
%                      ch_ref - vector contain valid channels (0 invalid, 1 valid)
%
%   bursts1, va, art will be the same length as the EEG data in ndm

% Read in EEG data from EDF file
%fname = '001_150411_000.edf';


fs2 = 64; %bval = 5; tl = 1; % Set sampling frequency to 64Hz, threshold for burst detector to 5 micovolts, and minimum burst duration to 1s
QQ = size(nref);
ch_left = nref(:,1)'; ch_right = nref(:,2)';

chr = [ch_left(ch_left~=0) ch_right(ch_right~=0)];
A = size(dat);
ndm = zeros(A); % Initialise matrix for raw EEG
reeg = zeros(A(1), floor((A(2)/(2*fs2)))); % Initialise matrix for rEEG
tst1 = zeros(1,A(1));
for z2 = 1:length(chr);
    ndm(chr(z2),:) = dat(chr(z2),:);
    reeg(chr(z2),:) =  estimate_rEEG(ndm(chr(z2),:), fs2);
    tst1(chr(z2)) = median(abs(hilbert(ndm(chr(z2),:))));
end
        
% SEARCH FOR ELECTRODE SHORTS
thr =  median(median(abs(hilbert(ndm(chr, :)))));
rf = find(tst1>thr/2);
%qq = chr(rf);
ch_left = rf(rf<=4);
ch_right = rf(rf>4);

% SEARCH FOR BURSTS   
bursts = cell(1, A(1)); 
for z3 = 1:length(chr);
%     [ba, ~, ~]=detector_per_channel_kirsi(ndm(chr(z3),:), fs2, tl, bval); % Using Kirsi's burst detector.
     ba=detector_per_channel_palmu(ndm(chr(z3),:), fs2, tst1(chr(z3))); % Using Kirsi's burst detector.

     bursts{1, chr(z3)} = ba;
end

% SEARCH FOR IMPEDENCE CHECKS
ic1 = zeros(1, length(ndm));
ic1(abs(ndm(chr(1),:))<0.01) = 1;
ic = process_ba_1ch(ic1, fs2, 1, 1, 1);
ic = ic(1:length(ndm));
clear ic1 

% SEARCH FOR HIGH AMPLITUDE ARTEFACT
art = find_artefacts_amp(ndm, bursts, ch_left, ch_right);
art = sum(art);
art(art<2) = 0;
art(art>=2) = 1;
art = art | ic;
%%%%%%%%%%%%%%%%%
%
% change for no art
%
%%%%%%%%%%%%%%%%%
%art = zeros(1,length(art));

% REFORMAT BURSTS FROM  CELL INTO MATRIX
bursts1 = zeros(size(ndm)); 
for z3 = 1:length(chr);
     bursts1(chr(z3),:) = bursts{chr(z3)};
end

