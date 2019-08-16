function fv1 = estimate_mse(data, fs, art)


fs1 = 128;
M = size(data);
[b,a] = butter(6, [1 20]/fs1); eeg = zeros(M(1), M(2)/(fs/fs1));
for ii = 1:M(1)
    eeg(ii,:) = resample(filtfilt(b, a, data(ii,:)), fs1, fs);
end


epl = 100*fs1;
block_no = floor(length(eeg)/epl); feat1 = zeros(block_no, 3);
scale = 1:20; val = zeros(block_no, 1);
for ii = 1:block_no
    r1 = (ii-1)*epl+1; r2 = r1+epl-1;
    r3  = (ii-1)*60*fs+1; r4=r3+60*fs-1;
    if sum(art(r3:r4))==0
    feats = zeros(M(1), 3);
    for z1 = 1:M(1)
        dat = eeg(z1, r1:r2);
        se = zeros(1, length(scale));
        for z2 = 1:length(scale)
            dum = conv(dat, ones(1, scale(z2)))/scale(z2);
            dum = dum(ceil(scale/2):scale(z2):end);
            se(z2) = SampEn(2, 0.2*std(dum), dum);
        end
        feats(z1, 1) = sum(se);
        feats(z1, 2) = max(se);
        feats(z1, 3) = sum(diff(se(1:5)));
    end
    val(ii) = 1;
    feat1(ii,:) = median(feats,1);
    end
end
fv1 = median(feat1(val==1,:), 1);

