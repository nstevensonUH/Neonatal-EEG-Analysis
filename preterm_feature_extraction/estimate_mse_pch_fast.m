function fvx = estimate_mse_pch_fast(data, fs, art)


fs1 = 128;
M = size(data);
[b,a] = butter(6, [1 20]/fs1); eeg = zeros(M(1), M(2)/(fs/fs1));
for ii = 1:M(1)
    eeg(ii,:) = resample(filtfilt(b, a, data(ii,:)), fs1, fs);
end
art = resample(double(art)', fs1, fs)';
art(art>0.5) = 1; art(art<=0.5) = 0;
epl = 100*fs1;
block_no = floor(length(eeg)/epl); 
scale = 1:20; 
%fv1 = zeros(1,M(1)) ; fv2 = fv1; fv3 = fv1;
%tic
%for z1 = 1:M(1)
    %z1
    feats = zeros(block_no, 3); val = zeros(1,block_no);
    for ii = 1:block_no
        r1 = (ii-1)*epl+1; r2 = r1+epl-1;
        if sum(art(r1:r2))==0
            dat = eeg(r1:r2);
            %se = zeros(1, length(scale));
            for z2 = 9 %1:length(scale)
                dum = conv(dat, ones(1, scale(z2)))/scale(z2);
                dum = dum(ceil(scale/2):scale(z2):end);
                se = SampEn(2, 0.2*std(dum), dum);
            end
            feats(ii, 1) = se;
            %feats(ii, 2) = max(se);
            %feats(ii, 3) = sum(diff(se(1:5)));
            val(ii) = 1;
        end
    end
    %if sum(val)>0.1*block_no
     %   fv1(z1) = median(feats(val==1,1));
        %fv2(z1) = median(feats(val==1,2));
        %fv3(z1) = median(feats(val==1,3));
    %else
    %    fv1(z1) = NaN; fv2(z1) = NaN; fv3(z1) = NaN;
    %end
%end
%toc
fvx = median(feats(val==1,1)); %[fv1 ; fv2 ; fv3];

