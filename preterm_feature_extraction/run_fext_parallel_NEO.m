function [out1, out2, out3] = run_fext_parallel_NEO(fname, filts) 
% filts = 1, run correction filters, filts = 0 do nothing

integerRange = 2^15-1;
NominalGain = 2365;
doubleRange  = 5e6/NominalGain/2;

%xx = strfind(fname, 'RAW');
%if isempty(xx)==1; cln_file = 0; else; cln_file = 1; end
%[data, fs1] = read_brm3_already_unzipped(fname);
[data, fs1] = read_brm3_no_scaling_v1(fname);
data = data .* (doubleRange / integerRange);
data = resample(data, 256, fs1);

if fs1>100
    [B1,A1] = butter(4,2*[42 58]/fs1,'stop'); % 50Hz notch
    % PREPROCESS - check for differences in machines - i.e. filters, artefacts
    % etc
    data = filter(B1, A1, data);
end

A = size(data);  
art = zeros(A); art1 = art; %art2 = art;
block_no = floor(A(1)/fs1);
for ii = 1:2
    for jj = 1:block_no
        r1 = (jj-1)*fs1+1; r2= r1+fs1-1;
        art1(r1:r2,ii) = max(abs(data(r1:r2,ii)));
%        art2(r1:r2,ii) = mean(abs(data(r1:r2,ii)));
    end            
end

art = art1 > 500 | art1 < 0.01; art = art';
for ii = 1:2
ra1 = find(diff([0 art(ii,:) 0])==1);
ra2 = find(diff([0 art(ii,:) 0])==-1);
ra1 = ra1-5*fs1; ra1(ra1<1)=1;
ra2 = ra2+5*fs1; ra2(ra2>length(art)) = length(art);
for jj = 1:length(ra1)
   art(ii,ra1(jj):ra2(jj))=1;
end
end

if filts == 1

Pxx = zeros(1,2^15+1);
data = data';
for ii = 1:2
    dta = data(ii,:);
    dta(art(ii,:)==1)=0;
    [Pdum, f] = pwelch(dta-mean(dta), hamming(2^16), 2^15, 2^16, fs1);
    Pxx = Pxx+Pdum';
end
fref = find(f>0.25 & f<4);
f_filt = f(fref(find(Pxx(fref) == max(Pxx(fref)))));

% Correct for different BRM versions.
filt_apply = 0;
if f_filt<0.66
    filt_apply = 1;
    [B2,A2] = butter(2,3/fs1, 'high'); % 1Hz HPF - new filter old filter was [B2,A2] = butter(4,2/fs1, 'high')
    data = filter(B2, A2, data')';
    A = size(data);  
    art = zeros(A); art1 = art; %art2 = art;
    block_no = floor(A(1)/fs1);
    for ii = 1:2
        for jj = 1:block_no
            r1 = (jj-1)*fs1+1; r2= r1+fs1-1;
            art1(r1:r2,ii) = max(abs(data(r1:r2,ii)));
        end            
    end

    art = art1 > 500 | art1 < 0.01; art = art';
    for ii = 1:2
    ra1 = find(diff([0 art(ii,:) 0])==1);
    ra2 = find(diff([0 art(ii,:) 0])==-1);
    ra1 = ra1-5*fs1; ra1(ra1<1)=1;
    ra2 = ra2+5*fs1; ra2(ra2>length(art)) = length(art);
    for jj = 1:length(ra1)
       art(ii,ra1(jj):ra2(jj))=1;
    end
    end
end

else
    
    data = data';

end

% EXTRACT FEATURES (1h epochs, 30 minute overlap)
fs2 = 64;
dat64 = resample(data', fs2, fs1)';
a64 = resample(double(art'), fs2, fs1)'; a64(a64>0.5)=1; a64(a64<1)=0;
durbins=[1 8 ; 8 16 ; 16 32 ; 32 64 ; 64 128 ; 128 256 ; 1 8192];
ob = [0.5 2 ; 2 4 ; 4 8 ; 8 13 ; 13 32];
dlim = 5; % CHANGE TO xlim burst fit, potential limit based on CDF analysis (35ms at the moment - half a fast spike)
epl1 = 3600*fs1; olap1 = epl1*0.5;
epl2 = 3600*fs2; olap2 = epl2*0.5;
block_no = floor(A(1)/olap1)-1; 
feat = cell(1,block_no); valid = zeros(2,block_no);
for ii = 1:block_no
    r1 = (ii-1)*olap1+1; r2 = r1+epl1-1;
    r3 = (ii-1)*olap2+1; r4 = r3+epl2-1;
    SC = zeros(1,2); amp1x = SC; bdm = SC; bdst = SC; fvx = SC; tp = SC; alpha = SC; fse90 = SC;
    amp1 = zeros(4,2); ibi = zeros(3,2); bds = ibi;
    B1 = zeros(2); B2 = B1; B3 = B1;
    sk = zeros(length(durbins), 2); kt = sk;
    ff = zeros(length(ob),2);
    for ch = 1:2
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
         else
             valid(ch, ii)=1;
       end
    end
    feat{ii} = [amp1x ; amp1 ; ibi ; bds ; sk ; kt ; B1 ; B2 ; B3 ; SC ; bdm ; bdst ; fvx ; log(tp) ; ff ; alpha ; fse90];     
end

out1 = feat;
out2 = [filts 0 sum((valid'))/(2*block_no)*100];
out3 = fname; 

   