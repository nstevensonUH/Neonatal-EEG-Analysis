function fv1 = estimate_path_length(data, fs, art)

M = size(data);
fs1 = 128;
[b,a] = butter(4, 32/(fs/2)); eeg = zeros(M(1), length(data)/(fs/fs1));

for ii = 1:M(1)
    eeg(ii,:) = resample(filter(b, a, data(ii,:)), fs1, fs);
end
fs1 = 128;
epl = 60*fs1;
block_no = floor(length(eeg)/epl); feat1 = zeros(block_no, 1);
for ii = 1:block_no
    r1 = (ii-1)*epl+1; r2 = r1+epl-1;
    r3  = (ii-1)*60*fs+1; r4=r3+60*fs-1;
    s1 = []; t1 = []; w1 = [];% s2 = []; t2 = []; w2 = [];
    if sum(art(r3:r4))==0
    for z1 = 1:M(1)
        for z2 = 1:M(1)
            if z1 > z2
                % ELEMENTS OF COHERENCE ESTIMATE PSD and CPSDs
            [Pxy, f] = cpsd(eeg(z1, r1:r2), eeg(z2, r1:r2), fs1, 96, fs, fs1);  
            fref = find(f>=4 & f<=8);
            [Pxx, ~] = pwelch(eeg(z1, r1:r2), fs1, 96, fs, fs1);  
            [Pyy, ~] = pwelch(eeg(z2, r1:r2), fs1, 96, fs, fs1);  
            % SURROGATE TESTING COMPONENT
%             zM = AAFTsur(eeg(z2, r1:r2),19);
%             As = 0;
%             for z3 = 1:19
%                 Ys = zM(:, z3)';
%                 [Pxys, ~] = cpsd(eeg(z1, r1:r2), Ys, 128, 96, 256, 128);  
%                 [Pyys, ~] = pwelch(Ys, 128, 96, 256, 128);  
%                 cohs = Pxys./(sqrt(Pxx.*Pyys));
%                 As = max([As max(abs(imag(cohs(fref))))]);
%             end
%           %  Take the maximum value of the imaginary component of the
            % coherence within a known frequency range
            %As = 0;
            coh = Pxy./(sqrt(Pxx.*Pyy)); fv = max(abs(imag(coh(fref))));
            % Test for coherence (no coherence no edge) and build bi-directional graph Matlab style
%             A2 = fv;
%             if fv>As 
%                 A1 = fv; 
%             else
%                 A1 = 0;
%             end
            s1 = [s1 z1]; t1 = [t1 z2]; w1 = [w1 fv];
%             s2 = [s2 z1]; t2 = [t2 z2]; w2 = [w2 A2];
            end
        end
    end
    
    B1 = graph(s1,t1,1./w1);
%    B2 = graph(s2,t2,1./w2);
    pl1 = zeros(length(s1),1); %pl2 = pl1;
    for jj = 1:length(s1)
        [~, pl1(jj)] = shortestpath(B1, s1(jj), t1(jj));
        %[~, pl2(jj)] = shortestpath(B2, s2(jj), t2(jj));
    end
    feat1(ii) = mean(pl1); %mean(pl(pl(:,1)>0,:));
    end
end
dum = median(feat1(feat1>0));
fv1 = dum(1);
%fv2 = dum(2);




