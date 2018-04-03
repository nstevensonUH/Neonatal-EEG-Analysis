function [dat1, fs, c1, c2, art1, art2] = read_data_montage_ema(filename, val);
%
%
%
%
%
%
%
%

[dat, ~, label, fs, scle, ~]  = read_edf(filename);

load filters_new;
fc = 0.5; Num = [1 -1]; Den =  [1 -(1- fc.*2.*pi./64)];
switch val
    case 1
        str = cell(8,2); 
        str{1,1} = 'Fp1'; str{1,2} = 'C3'; 
        str{2,1} = 'C3'; str{2,2} = 'O1'; 
        str{3,1} = 'Fp1'; str{3,2} = 'T3'; 
        str{4,1} = 'T3'; str{4,2} = 'O1'; 
        str{5,1} = 'Fp2'; str{5,2} = 'C4'; 
        str{6,1} = 'C4'; str{6,2} = 'O2'; 
        str{7,1} = 'Fp2'; str{7,2} = 'T4'; 
        str{8,1} = 'T4'; str{8,2} = 'O2'; 
        c1 = 1:4; c2 = 5:8; 
    case 2
       str = cell(4,2); 
        str{1,1} = 'Fp1'; str{1,2} = 'T3'; 
        str{2,1} = 'T3'; str{2,2} = 'O1'; 
        str{3,1} = 'Fp2'; str{3,2} = 'T4'; 
        str{4,1} = 'T4'; str{4,2} = 'O2'; 
        c1 = [1 2]; c2 = [3 4];
    case 3
         str = cell(4,2); 
        str{1,1} = 'Fp1'; str{1,2} = 'C3'; 
        str{2,1} = 'C3'; str{2,2} = 'O1'; 
        str{3,1} = 'Fp2'; str{3,2} = 'C4'; 
        str{4,1} = 'C4'; str{4,2} = 'O2';
        c1 = [1 2]; c2 = [3 4];  
    case 4
        str = cell(2,2); 
        str{1,1} = 'Fp1'; str{1,2} = 'T3'; 
        str{2,1} = 'Fp2'; str{2,2} = 'T4'; 
        c1 = 1; c2 = 2;
    case 5
        str = cell(2,2); 
        str{1,1} = 'T3'; str{1,2} = 'O1'; 
        str{2,1} = 'T4'; str{2,2} = 'O2'; 
        c1 = 1; c2 = 2;
 end
        
data_mont = cell(1,length(str));
art1 = zeros(length(str), 2); art2 = art1;
for jj = 1:length(str)
    ref1 = zeros(1,length(dat));
    ref2 = zeros(1,length(dat));
    for ii = 1:length(dat);
        ref1(ii) = length(findstr(label{ii}', str{jj,1})); 
        ref2(ii) = length(findstr(label{ii}', str{jj,2}));
    end
    
qq1 = find(ref1==1,1);
qq2 = find(ref2==1,1);
if length(dat{qq1})~=length(dat{qq2})
    data_mont{jj} = int16(zeros(1,length(dat{1})));
else
    art1(jj,1) = quantile(abs(hilbert(filter(Num, Den, resample(filter(Num_50, Den_50, double(dat{qq1})).*scle(qq1), 64, fs(qq1)) ))), 0.95);
    art1(jj,2) = quantile(abs(hilbert(filter(Num, Den, resample(filter(Num_50, Den_50, double(dat{qq2})).*scle(qq2), 64, fs(qq2)) ))), 0.95);
    art2(jj,1) = quantile(abs(hilbert(filter(Num, Den, resample(filter(Num_50, Den_50, double(dat{qq1})).*scle(qq1), 64, fs(qq1)) ))), 0.5);
    art2(jj,2) = quantile(abs(hilbert(filter(Num, Den, resample(filter(Num_50, Den_50, double(dat{qq2})).*scle(qq2), 64, fs(qq2)) ))), 0.5);

    data_mont{jj} = dat{qq1}-dat{qq2}; 
end
end

dat1 = zeros(length(data_mont),length(dat{1})/4);
for ii = 1:length(data_mont);
    dum = filter(Num_50, Den_50, double(data_mont{ii}).*scle(ii));
    dum = resample(dum, 64, fs(1));
    dat1(ii,:) = filter(Num, Den, dum);
end
fs = fs(1); 







