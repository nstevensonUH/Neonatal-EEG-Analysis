%

% load EEG
filename = 'demo.edf'; % not real EEG
[dat, hdr, label, fs, scle, offs] = read_edf(filename);

% Re-organize into montage (already done above)
% but is EDF is referential then something like this
% str = cell(8,2); 
% str{1,1} = 'Fp1'; str{1,2} = 'C3'; 
% str{2,1} = 'C3'; str{2,2} = 'O1'; 
% str{3,1} = 'Fp1'; str{3,2} = 'T3'; 
% str{4,1} = 'T3'; str{4,2} = 'O1'; 
% str{5,1} = 'Fp2'; str{5,2} = 'C4'; 
% str{6,1} = 'C4'; str{6,2} = 'O2'; 
% str{7,1} = 'Fp2'; str{7,2} = 'T4'; 
% str{8,1} = 'T4'; str{8,2} = 'O2'; 
% data_mont = cell(1,length(str));
% for jj = 1:length(str)
%     ref1 = zeros(1,length(dat));
%     ref2 = zeros(1,length(dat));
%     for ii = 1:length(dat)
%         ref1(ii) = length(findstr(label{ii}', str{jj,1})); 
%         ref2(ii) = length(findstr(label{ii}', str{jj,2}));
%     end
%     qq1 = find(ref1==1,1);
%     qq2 = find(ref2==1,1);
%     if length(dat{qq1})~=length(dat{qq2})
%         data_mont{jj} = int16(zeros(1,length(dat{1})));
%     else
%         data_mont{jj} = dat{qq1}-dat{qq2}; 
%     end
% end

fs1 = fs(1);
data_mont = dat; clear dat;
data = zeros(length(data_mont), length(data_mont{1}));
for ii = 1:length(data_mont)
    data(ii,:) = double(data_mont{ii}).*scle(ii);
end
% Filter (use whatever you want here)
[B1,A1] = butter(4,2*[42 58]/fs1,'stop'); % 50Hz notch filter
[B2,A2] = butter(4,2/fs1, 'high'); % maybe a high pass filter to eliminate slower oscillations
data = filtfilt(B1, A1, data')'; % filtfilt or filter can be used here, when in doubt use filter
data = filtfilt(B2, A2, data')';

% Detect Artefacts - in this case just look for bits that are way too low or way too high
% in this case 1 is an artefact and zero is no artefact if you want to skip
% this stage just use - A = size(data); art = zeros(A);
A = size(data);  
art = zeros(A); art1 = art; %art2 = art;
block_no = floor(A(2)/fs1);
for ii = 1:A(1)
    for jj = 1:block_no
        r1 = (jj-1)*fs1+1; r2= r1+fs1-1;
        art1(ii,r1:r2) = max(abs(data(ii,r1:r2)));
    end            
end
art = art1 > 500 | art1 < 0.01; 
for ii = 1:A(1)
    ra1 = find(diff([0 art(ii,:) 0])==1);
    ra2 = find(diff([0 art(ii,:) 0])==-1);
    ra1 = ra1-5*fs1; ra1(ra1<1)=1;
    ra2 = ra2+5*fs1; ra2(ra2>length(art)) = length(art);
    for jj = 1:length(ra1)
       art(ii,ra1(jj):ra2(jj))=1;
    end
end

% Extract Features
[feats, flist] = extract_features(data, art, fs1, 600);
