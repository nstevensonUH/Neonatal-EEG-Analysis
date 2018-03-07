function ema = estimate_ema(filename, val)
% This functions estimates the EEG maturational age from an EDF recording
% (filename).
%
% INPUTs: filename - the name of the EDF recording (ensure the EDF is saved as a referential recording with no filtering)
%                 val - a number that indexes different SVR models were generated for several different montages
% for variable, val
% 1 is an 8-channel montage (Fp-T, T-O, Fp-C, C-O)
% 2 is a 4-channel montage (Fp-T, T-O)
% 3 is a 4-channel montage (Fp-C, C-O)
% 4 is a 2-channel montage (Fp-T)
% 5 is a 2-channel montage (T-O)
%
% Aretfact detection based on the analysis of amplitude is also performed
% NANs or text will be output at run time stating the the EEG did not meet
% amplitude criteria or that the deisred montage could not be formulation
%
% Nathan Stevenson
% University of Helsinki, Finland
% January 2018

% Index of channels that are temporal for different SVR models
tht = cell(1,8); zz = tht;
tht{1} = [3 4 7 8];
tht{3} = [];
tht{2} = [1 2 3 4];
tht{4} = [1 2];
tht{5} = [1 2];

% Feature index for different SVR models
f1 = 1:46;
f2 = [1:15 17:38 40:46];
zz{1} = f1;
zz{3} = f2;
zz{2} = f1;
zz{4} = f1;
zz{5} = f1;

% READ FILE
try
    [dat, fs1, c1, c2, art1, art2] = read_data_montage_ema(filename, val);
catch
    disp('EEG montage could not be constructed')
    art1 = 1; art2 = 1;
end

if min(min(art1))<25 || min(min(art2)) < 0.5*median(median(art2)) 
    disp('EEG amplitude limits on at least 1 channel exceeded in referential montage')
    ema = NaN;
else

% PREPROCESS (detect bursts, sleep state and artefact) AND ESTIMATE FEATURES
ob = [0.5 3 ; 3 8 ; 8 15 ; 15 30]; % frequency bands of interest
len1 = length(dat)./fs1/3600; 
fs2 = 64; epl = 60; ep = epl*60*fs2; olap = ep/4; % segment into hour long epochs 15 minute overlap if possible
block_no = floor(length(dat)/olap)-3;
if block_no <= 0
    r1 = 1; r2 = length(dat);
    dd = dat(:, r1:r2);
    [bursts, art, ndm, reeg, ~, ~] = estimate_features_preprocessing(dd, fs1, [c1 ; c2]');    
    if sum(art)<12*60*fs2
ft1 = general_features_ss(ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
va = estimate_vigilance_state_new(bursts, c1, c2, 300*fs2, art);
ft2 = general_features_qs(va, ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
    else
        disp('Epoch 1 exceeds EEG amplitude limits')
        ft1=NaN*ones(1,23);
        ft2=NaN*ones(1,23);
    end
else
ft1 = zeros(block_no, 23); ft2 = ft1; %chr = zeros(1,block_no);
for z1 = 1:block_no;   
    r1 = 1+(z1-1)*olap; r2 = r1+ep-1;
    dd = dat(:, r1:r2);
    [bursts, art, ndm, reeg, ~, ~] = estimate_features_preprocessing(dd, fs1, [c1 ; c2]');    
     if sum(art)<12*60*fs2
    ft1(z1,:) = general_features_ss(ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
va = estimate_vigilance_state_new(bursts, c1, c2, 300*fs2, art);
ft2(z1,:) = general_features_qs(va, ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
     else
          disp(['Epoch ' num2str(z1) ' exceeds EEG amplitude limits'])
         ft1(z1,:)=NaN*ones(1,23);
        ft2(z1,:)=NaN*ones(1,23);
     end
end

end


% DO SVR
load(['model_file_' num2str(val)]);
ft = [ft1' ; ft2'];
ft = ft(zz{val},:);
A = size(ft); ema = zeros(1, A(2));
for ii = 1:A(2)
  ema(ii) = predict_svr(mdl1, ft(:,ii)', mu1, sig1);
end

end

