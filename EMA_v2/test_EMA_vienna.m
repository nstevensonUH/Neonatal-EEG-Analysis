% test space for developing standalon EMA estimator Vienna
%
%



addpath('E:\Neonatal_EEG_Data\preterm_Vienna\Matlab') 
addpath('E:\Neonatal_EEG_Data\preterm_Vienna') 
addpath('E:\Neonatal_EEG_Data\preterm_Vienna\new_EDFs_full') 
addpath('E:\Document_and_Code\papers\preterm_maturation\lib_svm\Matlab') 


% 
% 
% rr = [1 3 2 5 8];
% for ii = 3:5;
%     ii
% eval(['load test_run_' num2str(rr(ii))])
% %load test_run_1
% nref = 1:length(feat1);
% M = length(nref); stp = floor(M/3); nsp = [ 1:stp:M-2; stp:stp:M+stp-3]'; nsp(3,2)=M; 
% try
% [str1, perf] = optimise_svm_3f(ma1, feat1, [], 1:46, -2:4, -8:-4, [0.1 0.2 0.3 0.4 0.5], nsp);         
% [nx, nxt, nga, gat] = organise_data(ma1, feat1, nref, [], 1:46);
% [ny, mu, sig] = zscore(nx);
% mdl = svmtrain(nga, ny, str1);
% eval(['save model_file_' num2str(rr(ii)) ' mdl mu sig']);
% catch
% end
% end
% % B = size(nx); est_svm = zeros(B(1), 1);
% svl = mdl.totalSV;
% alpha = mdl.sv_coef;
% sv = mdl.SVs;
% gamma = mdl.Parameters(4);
% b = mdl.rho;
% for zn = 1:B(1);
%    val = (nx(zn,:)-mu)./sig;
%     [dd, ~, ~] = svmpredict(nga(zn), val, mdl, '-q');  
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     prd = zeros(1,svl);
%         for ii = 1:svl
%             prd(ii) = alpha(ii).*exp(-gamma.*norm(val-sv(ii,:)).^2);
%         end
%         %prd = prd-b;
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     est_svm(zn) = dd;
%     est_svm_sten(zn) = sum(prd)-b;
% end          
            
% NOTES
% REFS
% 1 - montage full 8 channel
% 2 - FpT and TO
% 3 - FpC and CO (no assessmen of temporal theta here)
% 4 - FpT
% 5 - TO
% linked to model file names.


filename = '001_250211_000.edf';
%SUBSEQUENT FEATURES
tht = cell(1,8); zz = tht;
tht{1} = [3 4 7 8];
tht{3} = [];
tht{2} = [1 2 3 4];
tht{4} = [1 2];
tht{5} = [1 2];

f1 = 1:46;
f2 = [1:15 17:38 40:46];
zz{1} = f1;
zz{3} = f2;
zz{2} = f1;
zz{4} = f1;
zz{5} = f1;

% READ FILE
prd1= zeros(5,2);
for qq = 1:5
val = qq;
filename = '001_250211_000.edf';
try
    [dat, fs1, c1, c2, art1] = read_data_montage_ema(filename, val);
catch
    disp('EEG montage could not be constructed')
end

if art1


% PREPROCESS AND ESTIMATE FEATURES (change current filename)
ob = [0.5 3 ; 3 8 ; 8 15 ; 15 30];
len1 = length(dat)./fs1/3600; 
fs2 = 64; epl = 60; ep = epl*60*fs2; olap = ep/4;
block_no = floor(length(dat)/olap)-3;
if block_no == 0
    r1 = 1; r2 = length(dat);
    dd = dat(:, r1:r2);
    [bursts, art, ndm, reeg, ch_left, ch_right] = estimate_features_preprocessing(dd, fs1, [c1 ; c2]');    
ft1 = general_features_ss(ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
va = estimate_vigilance_state_new(bursts, c1, c2, 300*fs2, art);
ft2 = general_features_qs(va, ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});

else
ft1 = zeros(block_no, 23); ft2 = ft1; %chr = zeros(1,block_no);
for z1 = 1:block_no;   
    r1 = 1+(z1-1)*olap; r2 = r1+ep-1;
    dd = dat(:, r1:r2);
    [bursts, art, ndm, reeg, ch_left, ch_right] = estimate_features_preprocessing(dd, fs1, [c1 ; c2]');    
    ft1(z1,:) = general_features_ss(ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
va = estimate_vigilance_state_new(bursts, c1, c2, 300*fs2, art);
ft2(z1,:) = general_features_qs(va, ndm, bursts, reeg, art, c1, c2, ob, fs2, [c1' c2'], tht{val});
end

end


% DO SVR
load(['model_file_' num2str(val)]);
ft = [ft1' ; ft2'];
ft = ft(zz{val},:);
A = size(ft);
for ii = 1:A(2)
  prd1(qq, ii) = predict_svr(mdl, ft(:,ii)', mu, sig);
end

end

filename = '001_250211_000.edf';
for qq = 1:5
    ema = estimate_ema(filename, qq);
end

