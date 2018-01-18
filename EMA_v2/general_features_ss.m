function  fv1 = general_features_ss(ndm, bursts, reeg, art, ch_left, ch_right, ob, fs2, prs, tht);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fv1 = zeros(1, 23);
ch = [ch_left ch_right];
% Estimate amplitude envelope
val1 = zeros(1, length(ch)); val2 = val1; val3 = val1;
for z1 = 1:length(ch)
     aa = hilbert(ndm(ch(z1),:));
     vv = quantile(abs(aa(art==0)), [0.5 0.05 0.95]); % Ignore artefact, estimate summary stats per channel
     val1(z1) = vv(1); val2(z1) = vv(2); val3(z1) = vv(3);
end

% Summary stats across channels
fv1(1) = median(val1);
fv1(2) = median(val2); 
fv1(3) = median(val3); 

% Estimate rEEG
a1 = find(diff([0 art 0])==1); a2 = find(diff([0 art 0])==-1);
a1 = ceil(a1/128); a2 = floor(a2/128); % Resample artefact annotation
aref = zeros(1,length(reeg)); 
for jj = 1:length(a1); aref(a1(jj):a2(jj))=1; end
rr = reeg(ch, aref==0); % Ignore artefact
vv = quantile(rr', [0.5 0.05 0.95]); % summarise per channel
fv1(4) = median(vv(1,:));  % summary stats across channels
fv1(5) = median(vv(2,:)); 
fv1(6) = median(vv(3,:)); 

clear val rr 

% Estimate features based on the burst/SAT annotation
N1 = length(ndm(1,:));
a1 = find(diff([0 1-art 0]) == 1); a2 = find(diff([0 1-art 0]) == -1); if a2(end)==N1+1; a2(end)=N1; end
bb1 = zeros(1, length(ch)); bb2 = bb1; bb3 = bb1; bb4 = bb1; bb5 = bb1;
ib1 = zeros(1, length(ch)); ib2 = ib1; ib3 = ib1; ib4 = ib1;
for z1 = 1:length(ch);
    ba = bursts(ch(z1),:); bd = []; ibd = [];
    for z2 = 1:length(a1)
        r1 = find(diff([0 ba(a1(z2):a2(z2)) 0])==1);
        r2 = find(diff([0 ba(a1(z2):a2(z2)) 0])==-1)-1;
        ibd = [ibd r1(2:end)-r2(1:end-1)];
        if ba(a1(z2))==1; r1 = r1(2:end); r2 = r2(2:end); end
        if ba(a2(z2))==1; r1 = r1(1:end-1); r2 = r2(1:end-1); end
        bd = [bd r2-r1];
%         for z3 = 1:length(r1);
%             b_area = [b_area mean(abs(hilbert(ndm(ch(z1), r1(z3):r2(z3))))).*((r2(z3)-r1(z3))/fs2)];
%         end 
    end  
    bb = quantile(bd, [0.5 0.05 0.95]);   % summarise burst duration per channel
    bb1(z1) = bb(1)/fs2; bb2(z1) = bb(2)/fs2; bb3(z1) = bb(3)/fs2; 
    bb4(z1) = sqrt(mean((bd./fs2).^2)); % RMS burst duration
    bb5(z1) = length(bd)/(sum(art==0)/fs2/3600);  % Burst number per hour
    ib = quantile(ibd, [0.5 0.05 0.95]); % summarise inter-burst duraton per channel
    ib1(z1) = ib(1)/fs2; ib2(z1) = ib(2)/fs2; ib3(z1) = ib(3)/fs2; 
    ib4(z1) = sqrt(mean((ibd./fs2).^2)); % RMS inter-burst duration
    
end

% Summarise across channels
fv1(7) = median(bb5(isnan(bb1)==0)); 
fv1(8) = median(bb4(isnan(bb1)==0));
fv1(9) = median(bb1(isnan(bb1)==0));
fv1(10) = median(bb2(isnan(bb1)==0));
fv1(11) = median(bb3(isnan(bb1)==0));
fv1(12) = median(ib4(isnan(ib1)==0));
fv1(13) = median(ib1(isnan(ib1)==0)); 
fv1(14) = median(ib2(isnan(ib1)==0)); 
fv1(15) = median(ib3(isnan(ib1)==0)); 

% Temporal theta
f = 0:fs2/N1:fs2-1/N1;
DTA = zeros(1,N1);
%qq = find(ch == 3  | ch == 4 | ch == 7  | ch == 8); % Must select temporal channels
qq = [];
for z1=1:length(tht);
    qrf = find(ch==tht(z1));
    qq = [qq qrf];
end

for z1 = 1:length(qq)
     dta = ndm(ch(qq(z1)), :);
     dta(art==1) = 0;
     dta = dta-mean(dta);
     N = sum(art==0);
     DTA = DTA+1/N*abs(fft(dta)).^2;
end
DTA = 2*DTA(1:ceil(N1/2))/length(ch);
fref1 = find(f>=3 & f<8);
fv1(16) = mean(DTA(fref1)); 

% IGNORE SYNCHRONY MEASURES FOR NOW
load iir_filter_05_16_64
% % ASI and other Synchrony
%prs = [1 5 ; 2 6 ; 3 7 ; 4 8];
D = size(prs);
c1 = 1; % These pairs will change.
r_p_env = []; asi = [];
for ii = 1:D(1);
    p1 = prs(ii,:);
    if isempty(find(ch_left == p1(1)))==0 & isempty(find(ch_right == p1(2)))==0
        asi(c1) = asi_estimate(ndm(p1(1),:), ndm(p1(2), :), art, fs2);
%        r_sp(c1) = corr(bursts{1, p1(1)}', bursts{1, p1(2)}', 'type', 'Spearman');
      %   slike(c1) = SynchonizationLikelihood(ndm(p1(1),:), ndm(p1(2), :), fs2, 0.5, 16);
        h1 = abs(hilbert(filter(Numl, Denl, filter(Numh, Denh, ndm(p1(1),:)))));
         h2 = abs(hilbert(filter(Numl, Denl, filter(Numh, Denh, ndm(p1(2),:)))));
         r_p_env(c1) = corr((h1(art==0))', (h2(art==0))', 'type', 'Pearson');
         c1 = c1+1;
     end
end
fv1(17) = median(asi);
fv1(18) = median(r_p_env); 

% SPECTRAL ANALYSIS
f = 0:fs2/N1:fs2-1/N1;
DTA = zeros(length(ch),N1);
BB = size(ob);
% Estimate PSD - periodogram across all channels
for z1 = 1:length(ch)
     dta = ndm(ch(z1), :);
     dta(art==1) = 0;
     dta = dta-mean(dta);
     N = sum(art==0);
     DTA(z1,:) = 1/N*abs(fft(dta)).^2; % periodogram (slighlty modified due to artefact)
end
DTA = 2*DTA(:,1:ceil(N1/2)); % still works as DTA(1) will be zero

fref1 = find(f>=0.5 & f<30); % full band of interest
tp = sum(DTA(:,fref1)');
fv1(19) = log(median(tp)); % ~signal energy

% Estimate relative spectral powers in each band
ff = zeros(length(ch), BB(1));
for z1 = 1:length(ch)
    for z2 = 1:BB(1)
        fref1 = find(f>=ob(z2,1) & f<ob(z2, 2)); 
        ff(z1, z2) = sum(DTA(z1, fref1))./tp(z1);
    end
end
% Summarise across channels
fv1(20) = median(ff(:,1));
fv1(21) = median(ff(:,2));
fv1(22) = median(ff(:,3));
fv1(23) = median(ff(:,4));
end