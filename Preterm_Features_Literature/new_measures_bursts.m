function fv = new_measures_bursts(dat, prs, fs, art)
% Estimate several features from the literature that show correlation with
% PMA in preterm infants (<38 weeks GA)
%
% dat - EEG data (M x N), M - number of recording channels and N length of
% recording in samples
% prs - index defining which channels are opposed to each other (or which
% channels you wish to compare to each other), e.g. prs = [1 2 ; 3 4] -
% will perform cross-channel comparisons, in the dat variable, between 1
% and 2, and 3 and 4, and then everage the results.
% fs - sampling frequency
% art - artefact mask, 1 for artefact zero for no-artefact, if no artefact
% detectio is applied then let art = zeros(1, N).
%
% Nathan Stevenson
% QIMR Berghofer 2019
%

chr = sort(prs(:));
data = dat;
epl = fs*2; N = length(data(1,:));
block_no = floor(N/epl); art_reeg = zeros(length(chr), block_no);
for ii = 1:block_no
    q1 = (ii-1)*epl+1; q2 = q1+epl-1;
    for kk = 1:length(chr)
      art_reeg(chr(kk),ii) = max(art(q1:q2));
    end
end
fs2 = 100;
[B, A] = butter(4, 1/fs2, 'high');
M = 2*fs2; 
d0 = 0.15;
f = 0:fs/M:fs-1/M;
fref = find(f>=0 & f<=13.5);
dat1 = resample(dat', fs2, fs)';
art2 = resample(double(art), fs2, fs); art2(art2>0.5)=1; art2(art2<=0.5)=0;
dat1(:, art2==1) = NaN;
tref = 200:M:length(dat)-1;

% Segment Rate
SR = zeros(1,length(chr));
for kk = 1:length(chr)
    dat1 = filter(B,A, dat1);
    block_no = floor(length(dat1)/M);
    D = zeros(1, block_no-1); dval = zeros(1, block_no-1);
    for ii = 1:block_no-1
        r1 = (ii-1)*M+1; r2 = r1+M-1;
        r3 = ii*M+1; r4 = r3+M-1;
        if sum(isnan(dat1(r1:r2)))>1 || sum(isnan(dat1(r3:r4)))>1
            dval(ii) = 1;
        else
        X1 = abs(fft(dat1(r1:r2)));
        X2 = abs(fft(dat1(r3:r4)));
        [~,~,D(ii)] = kstest2(X1(fref), X2(fref));
        end
    end
    sref = find(D>d0*2); % can select to maximise variability
    SR(kk) = length(sref)/sum(dval==0);
end

fv(1) = median(SR);

% rEEG
for kk = 1:length(chr)
    reeg =  estimate_rEEG(data(chr(kk),:), fs);
    bw(kk) = diff(quantile(reeg(isnan(art_reeg(chr(kk),:))==0), [0.05 0.95]));
end
fv(2) = median(bw);

% Suppression Curve
data = dat;
[B, A] = butter(4, 1/fs, 'high');
data = filter(B,A, data')';
fv(3)=calculateSC_NS(data(chr,:), 256, art);

% Global and hemispheric ASI
[fv(4), fv(5)] = estimate_ASI(data, fs, prs, art);

% mean phase locking index
fv(6) = estimate_mPLI(data(chr,:), fs, art);
% path length (coherence)
fv(7) = estimate_path_length(data(chr,:), fs, art);
% multi-scale entropy (slope first 5 scales, mean and maxiumu)
fvx = estimate_mse(data(chr,:), fs, art);
fv(8) = fvx(1);
fv(9) = fvx(2);
fv(10) = fvx(3);

