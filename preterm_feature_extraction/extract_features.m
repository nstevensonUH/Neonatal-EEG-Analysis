function [feats, flist] = extract_features(data, art, fs1, epl)
% This function extracts a series of features from a period of filtered EEG
%
% Inputs: data is the filtered EEG epoch you wish to analyse, it should be
%              of size channel x time, e.g. an 8 channel recording for 1 minute at 250Hz
%              will mean size(data) = [8,15000]
%         art is an artefact detection matrix
%         fs1 is the sampling frequency of the input data
%         epl is the epoch/segment length (in seconds) over which you want to estiamte
%             the features, a 50% overlap is used if the data is longer than
%             epl (note for burst features long epochs are required, i.e. 1
%             hour, there is also a lower limit here, say nothing less than
%             1 minute
%           
%
% Outputs: feats is a cell array that is size [1, M], where M is the number of segments 
%                of length epl that are extracted from the input data
%                matrix. Each cell entry contains a matrix of size [43, K],
%                where 43 is the number of features estimated per
%                epoch/segment and K is the number of EEG channels.
%          flist is a cell that is size [1,43] containing a text descriptor
%                of each feature (for large data process it is probably best to ingore this output)
%        
% Dependencies: calculateSC_NS_pch (this calculates the suppression curve see Dereymaeker A, Koolen N, 
% Jansen K, Vervisch J, Ortibus E, De Vos M, Van Huffel S, Naulaers G. The suppression curve as a 
% quantitative approach for measuring brain maturation in preterm infants. Clinical Neurophysiology. 
% 2016 Aug 1;127(8):2760-5.),
% detector_per_channel_palmu (a burst detector see Palmu K, Stevenson N, Wikström S, Hellström-Westas L,
% Vanhatalo S, Palva JM. Optimization of an NLEO-based algorithm for automated detection of spontaneous
% activity transients in early preterm EEG. Physiological measurement. 2010 Oct 11;31(11):N85),
% estimate_mse_pch_fast (the multi-scale entropy calculates sample entropy over a range
% of subsampled versions of the original signal, this is excessive and we just calculate 
% the sample entropy for the optimal subsampling for correlating with age in preterm infants)
% It also the Matlab rmoutliers function which is only available on newer
% versions of Matlab from 2018b
% 
% Nathan Stevenson
% QIMR Berghofer
% May 2022

% EXTRACT FEATURES (1h epochs, 30 minute overlap)
AA = size(data);
fs2 = 64;
dat64 = resample(data', fs2, fs1)'; % Almost everything is calculated on this except Sample Entropy
a64 = resample(double(art'), fs2, fs1)'; a64(a64>0.5)=1; a64(a64<1)=0;
durbins=[1 8 ; 8 16 ; 16 32 ; 32 64 ; 64 128 ; 128 256 ; 1 8192]; % durations in samples
ob = [0.5 2 ; 2 4 ; 4 8 ; 8 13 ; 13 32]; % frequency bands in Hz (we split delta into two subbands)
dlim = 4; % CHANGE TO xlim burst fit, potential limit based on CDF analysis (70ms at the moment - half a fast spike)
epl1 = epl*fs1; olap1 = epl1*0.5;
epl2 = epl*fs2; olap2 = epl2*0.5;
block_no = floor(AA(2)/olap1)-1; 
feats = cell(1,block_no); 
for ii = 1:block_no % loop per block
    r1 = (ii-1)*olap1+1; r2 = r1+epl1-1;
    r3 = (ii-1)*olap2+1; r4 = r3+epl2-1;
    SC = zeros(1,AA(1)); amp1x = SC; bdm = SC; bdst = SC; fvx = SC; tp = SC; alpha = SC; fse90 = SC;
    amp1 = zeros(4,AA(1)); ibi = zeros(3,AA(1)); bds = ibi;
    B1 = zeros(2); B2 = B1; B3 = B1;
    sk = zeros(length(durbins), 2); kt = sk;
    ff = zeros(length(ob),2);
    for ch = 1:AA(1) % loop per channel
         dat = dat64(ch, r3:r4);
         dat1 = data(ch, r1:r2);
         adat = a64(ch, r3:r4);
         adat1 = art(ch, r1:r2);
         if sum(adat)<0.5*length(adat)
         % AMPLITUDE
         dum = quantile(abs(hilbert(dat)), [0.05 0.25 0.5 0.75 0.95]);
         amp1x(ch) = dum(3);
         amp1(:,ch) = dum([1 2 4 5])';
         % BURST DURATION & INTERVAL
         ba = detector_per_channel_palmu(dat, fs2, amp1x(ch)); % Using Kirsi's burst detector. % CHECK FOR SPEED
         r1z = find(diff([0 ba 0]) == 1);
         r2z = find(diff([0 ba 0]) == -1);
         ibis = r1z(2:end)-r1z(1:end-1);
         bdurs = r2z-r1z;

         a1z = find(diff([0 adat 0]) == 1);
         a2z = find(diff([0 adat 0]) == -1);
         refx = 1:length(bdurs);
         for pp = 1:length(a1z)
             ref1 = [find(abs(r1z-a1z(pp)) == min(abs(r1z-a1z(pp))),1) find(abs(r2z-a1z(pp)) == min(abs(r2z-a1z(pp))),1)];
             ref2 = [find(abs(r1z-a2z(pp)) == min(abs(r1z-a2z(pp))),1) find(abs(r2z-a2z(pp)) == min(abs(r2z-a2z(pp))),1)];
             val = [min([ref1 ref2])-1 max([ref1 ref2])+1];    
             if val(1)<1; val(1)=1; end
             if val(end)>length(bdurs); val(end)=length(bdurs); end
             refx(val(1):val(2))=0;
         end
         ibi(:,ch) = quantile(ibis(refx(1:end-1)>0), [0.05 0.5 0.95])/fs2;
         bds(:,ch) = quantile(bdurs(refx>0), [0.05 0.5 0.95])/fs2;

         % SKW / KURT BURSTS 
        amp = abs(hilbert(dat)).^2;
        amp(adat==1)=0;
        amp = conv(amp, [1 1 1 1 1]/5, 'same');
        th1 = quantile(amp, 100); r1q = zeros(1,length(th1));
        for zz = 1:length(th1)
            dum = zeros(1, length(amp));
            dum(amp<th1(zz)) = 0;
            dum(amp>=th1(zz)) = 1;      
            r1x = find(diff([0 dum 0]) == 1);
            r2x = find(diff([0 dum 0]) == -1);
            tst2 = r2x-r1x;
            rf = find(tst2<=dlim);
            for z2 = 1:length(rf)
                dum(r1x(rf(z2)):r2x(rf(z2)))=0;
            end
            r1q(zz) = length(find(diff([0 dum 0])==1));
        end
        th = th1(find(r1q==max(r1q),1)); 
        dum(amp<th) = 0;
        dum(amp>=th) = 1;
        r1x = find(diff([0 dum 0]) == 1);
        r2x = find(diff([0 dum 0]) == -1);
        tst2 = r2x-r1x;
        r1x = r1x(tst2>dlim);
        r2x = r2x(tst2>dlim);
        tst2 = tst2(tst2>dlim);
        M = ceil(2*mean(tst2)*10); 
        for qq = 1:length(durbins)
            bav = zeros(1,M);
            r1y = r1x(tst2 > durbins(qq,1) & tst2 <= durbins(qq,2));
            r2y = r2x(tst2 > durbins(qq,1) & tst2 <= durbins(qq,2));
            for kk = 1:length(r1y)
               dum =  amp(r1y(kk):r2y(kk)-1)-th; 
               xx = linspace(1,length(dum),M);
               pp = pchip(1:length(dum), dum, xx); 
               pp = pp./sum(pp);
               bav = bav+pp; % Average shape
            end
            % 4 moments of burst shape only need skewness and kurtosis
            bav1 = (bav-min(bav)); bav1 = bav1./sum(bav1);
            xx = linspace(0,1,M);
            mn = sum(xx.*bav1);
            sd = sqrt(sum((xx-mn).^2.*bav1));
            sk(qq,ch) = sum((xx-mn).^3.*bav1)./sd^3;
            kt(qq,ch) = sum((xx-mn).^4.*bav1)./sd^4-3;
        end

        bd1 = zeros(1,length(r1x)); ba1 = bd1;
        for kk= 1:length(r1x)
            bd1(kk) = length(r1x(kk):r2x(kk)-1)/fs2;
            ba1(kk) = trapz(amp(r1x(kk):r2x(kk)-1)-th)/fs2; % Area Trapezoidal
        end
        B1(:,ch) = polyfit(log(bd1), log(ba1), 1);
        B2(:,ch) = polyfit(durbins(:,1)'./fs2, sk(:,ch)', 1);
        B3(:,ch) = polyfit(durbins(:,1)'./fs2, kt(:,ch)', 1);      
        SC(ch) = calculateSC_NS_pch(dat, fs2, adat);
        bdm(ch) = mean(tst2)/fs2;
        bdst(ch) = std(tst2)/fs2;

        % Sample Entropy
        fvx(ch) = estimate_mse_pch_fast(dat1, fs1, adat1);

        % SPECTRAL ANALYSIS, POWER/SLOPE/SPECTRAL EDGE (90%)
        N1 = length(dat);
        f = 0:fs2/N1:fs2-1/N1;
        BB = size(ob);
        % Estimate PSD - periodogram across all channels
        dta = dat; % use 64Hz sampled data
        dta(adat==1) = 0;
        dta = dta-mean(dta);
        N = sum(adat==0);
        DTA = 1/N*abs(fft(dta)).^2; % periodogram (slighlty modified due to artefact)
        DTA = 2*DTA(:,1:ceil(N1/2)); % still works as DTA(1) will be zero
        fref1 = find(f>=min(min(ob)) & f<max(max(ob))); % full band of interest
        tp(ch) = sum(DTA(fref1));
        % Estimate relative spectral powers in each band
        %ff = zeros(1,BB(1));
        for z2 = 1:BB(1)
           fref1 = find(f>=ob(z2,1) & f<ob(z2, 2)); 
           ff(z2, ch) = sum(DTA(fref1))./tp(ch);
        end
        tp(ch) = log(tp(ch));
        vv = floor(log2(length(dta)));
        [Pxx, f] = pwelch(dta, hamming(2^(vv-3)), 2^(vv-4), 2^vv-3, fs2);
        fr = 2.^[log2(min(min(ob))+1.5):0.1:log2(16)]; P = zeros(1,length(fr)-1);
        for cc = 1:length(fr)-1
           rx = find(f>=fr(cc) & f<(fr(cc+1)));
           P(cc) = mean(Pxx(rx)); 
        end  
       outliers = 1; logf = log2(fr); logP = log2(P); idx = zeros(1,length(P));
       while outliers ~= 0
            x1 = logf(idx==0); y1 = logP(idx==0);
            B = regress(y1', [x1' ones(size(x1'))]);
            res = y1 - (B(1)*x1+B(2)*ones(1,length(x1)));
            [~, idx] = rmoutliers(res);
            outliers = sum(idx);
       end
       alpha(ch) = B(1);
       fdum = f(f>=0.5 & f<=32);
       frx = find(cumsum(Pxx(f>=0.5 & f<=32))./sum(Pxx(f>=0.5 & f<=32))>0.9, 1); 
       fse90(ch) = fdum(frx); 
       end
    end
    feats{ii} = [amp1x ; amp1 ; ibi ; bds ; sk ; kt ; B1 ; B2 ; B3 ; SC ; bdm ; bdst ; fvx ; tp ; ff ; alpha ; fse90];     
end

flist{1} = 'Median Amplitude Envelope'; %Amplitude Envelope is estimate using the analytic associate of a signal (via Hilbert transform)
flist{2} = '5th Percentile Amplitude Envelope';
flist{3} = '25th Percentile Amplitude Envelope';
flist{4} = '75th Percentile Amplitude Envelope';
flist{5} = '95th Percentile Amplitude Envelope';
flist{6} = '5th Percentile Inter-burst Interval';
flist{7} = '50th Percentile Inter-burst Interval';
flist{8} = '95th Percentile Inter-burst Interval';
flist{9} = '5th Percentile Burst Duration';
flist{10} = '50th Percentile Burst Duration';
flist{11} = '95th Percentile Burst Duration';
flist{12} = 'Burst Skewness/Symmetry (15.625-125ms)';
flist{13} = 'Burst Skewness/Symmetry (125-250ms)';
flist{14} = 'Burst Skewness/Symmetry (250-500ms)';
flist{15} = 'Burst Skewness/Symmetry (0.5-1s)';
flist{16} = 'Burst Skewness/Symmetry (1-2s)';
flist{17} = 'Burst Skewness/Symmetry (2-4s)';
flist{18} = 'Burst Skewness/Symmetry (15.625ms-128s)';
flist{19} = 'Burst Kurtosis/Sharpness (15.625-125ms)';
flist{20} = 'Burst Kurtosis/Sharpness (125-250ms)';
flist{21} = 'Burst Kurtosis/Sharpness (250-500ms)';
flist{22} = 'Burst Kurtosis/Sharpness (0.5-1s)';
flist{23} = 'Burst Kurtosis/Sharpness (1-2s)';
flist{24} = 'Burst Kurtosis/Sharpness (2-4s)';
flist{25} = 'Burst Kurtosis/Sharpness (15.625ms-128s)';
flist{26} = 'Slope Burst Duration vs Burst Area';
flist{27} = 'Intercept Burst Duration vs Burst Area';
flist{28} = 'Slope Burst Duration vs Skewness';
flist{29} = 'Intercept Burst Duration vs Skewness';
flist{30} = 'Slope Burst Duration vs Kurtosis';
flist{31} = 'Intercept Burst Duration vs Kurtosis';
flist{32} = 'Suppression Curve';
flist{33} = 'Mean Burst Duration';
flist{34} = 'Standard Deviation Burst Durations';
flist{35} = 'Sample Entropy';
flist{36} = 'Log PSD Energy (wideband)';
flist{37} = 'Relative Delta 1 power';
flist{38} = 'Relative Delta 2 power';
flist{39} = 'Relative Theta power';
flist{40} = 'Relative Alpha power';
flist{41} = 'Relative Beta power';
flist{42} = 'Spectral Slope';
flist{43} = 'Spectral Edge Frequency (90th Percentile)';
   